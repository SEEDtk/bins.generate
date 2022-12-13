/**
 *
 */
package org.theseed.bins;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.p3api.P3Connection;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.Jsoner;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object manages a set of bins.  The bins are kept in a master list and an index is kept of which
 * bin contains each contig.
 *
 * @author Bruce Parrello
 *
 */
public class BinGroup implements Iterable<Bin> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinGroup.class);
    /** map of contig IDs to bins */
    private Map<String, Bin> contigMap;
    /** set of bins */
    private SortedSet<Bin> binList;
    /** statistics about the bin group */
    private CountMap<String> stats;
    /** contig file name */
    private File inputFile;
    /** default value for count map */
    private static final JsonObject noCounts = new JsonObject();
    /** default value for empty bin list */
    private static final JsonArray noBins = new JsonArray();

    private static enum GroupKeys implements JsonKey {
        COUNTS(noCounts),
        IN_FILE(null),
        BINS(noBins);

        /** default value for key */
        private Object mValue;

        private GroupKeys(Object value) {
            this.mValue = value;
        }

        @Override
        public String getKey() {
            return this.name().toLowerCase();
        }

        @Override
        public Object getValue() {
            return this.mValue;
        }

    }

    /**
     * Create a new, empty bin group.
     */
    public BinGroup() {
        this.setup(1000);
    }

    /**
     * Initialize the bin group with a contig hash the specified size.
     *
     * @param hashSize		expected number of contigs
     */
    protected void setup(int hashSize) {
        this.contigMap = new HashMap<String, Bin>(hashSize * 4 / 3 + 1);
        this.binList = new TreeSet<Bin>();
        this.stats = new CountMap<String>();
        this.inputFile = null;
    }

    /**
     * Load a bin group from a JSON file.
     *
     * @param inFile	file containing the bin group in JSON format
     *
     * @throws IOException
     * @throws JsonException
     */
    public BinGroup(File inFile) throws FileNotFoundException, IOException, JsonException {
        log.info("Loading bin group from {}.", inFile);
        try (Reader fileReader = new FileReader(inFile)) {
            // Read in the bin group.
            JsonObject groupObject = (JsonObject) Jsoner.deserialize(fileReader);
            // Get the bin list.
            JsonArray binArray = (JsonArray) groupObject.getCollectionOrDefault(GroupKeys.BINS);
            // Initialize this object's structures to hold the expected number of contigs.
            this.setup(binArray.size() * 2);
            // Add all the bins found.
            for (Object binObject : binArray) {
                Bin newBin = new Bin((JsonObject) binObject);
                this.addBin(newBin);
            }
            // Get the input file name.
            String fileString = groupObject.getStringOrDefault(GroupKeys.IN_FILE);
            if (fileString == null)
                this.inputFile = null;
            else
                this.inputFile = new File(fileString);
            // Read in the counts.
            JsonObject counts = groupObject.getMapOrDefault(GroupKeys.COUNTS);
            for (var countName : counts.keySet()) {
                int countValue = P3Connection.getInt(counts, countName);
                this.stats.count(countName, countValue);
            }
        }
        log.info("{} contigs and {} bins read from {}.", this.contigMap.size(), this.binList.size(), inFile);
    }

    /**
     * Load a bin group from a FASTA file.  We will compute the coverage here, and output the contigs eligible for the
     * seed-protein search to the specified output file.
     *
     * @param fastaFile		input file containing contigs from an assembly
     * @param parms			tuning parameters for contig filtering
     * @param reducedFile	output file for seed-search contigs
     *
     * @throws IOException
     */
    public BinGroup(File fastaFile, BinParms parms, File reducedFile) throws IOException {
        // Initialize the object structures.
        this.setup(1000);
        // Create the coverage filter.
        ContigFilter filter = new ContigFilter(parms);
        // Open the output file and connect to the input file.
        try (FastaInputStream inStream = new FastaInputStream(fastaFile);
                FastaOutputStream outStream = new FastaOutputStream(reducedFile)) {
            this.inputFile = fastaFile;
            log.info("Reading contigs from {}.", fastaFile);
            int seedUsableCount = 0;
            for (Sequence seq : inStream) {
                Bin seqBin = filter.computeBin(seq, this.stats);
                if (seqBin.getStatus() == Bin.Status.SEED_USABLE) {
                    // Here the sequence is good enough for the seed-protein search.
                    outStream.write(seq);
                    seedUsableCount++;
                }
                if (seqBin.getStatus() != Bin.Status.BAD) {
                    // Here the sequence is good enough to add to the bin group.
                    this.addBin(seqBin);
                }
            }
            log.info("{} seed-search sequences written to {}, {} saved for binning.",
                    seedUsableCount, this.binList.size());
        }
    }

    /**
     * Add a bin to this group.
     *
     * @param bin		bin to add
     */
    public void addBin(Bin bin) {
        this.binList.add(bin);
        for (String member : bin.getContigs())
            this.contigMap.put(member, bin);
    }

    /**
     * Merge two bins.
     *
     * @param bin		target bin
     * @param bin2		bin to merge into the first bin
     */
    public void merge(Bin bin1, Bin bin2) {
        // Denote bin2 is no longer in the group.
        this.binList.remove(bin2);
        // Combine the contigs.
        bin1.merge(bin2);
        // Point all the second bin's contigs to the first bin.
        for (String member : bin2.getContigs())
            this.contigMap.put(member, bin1);
    }

    /**
     * A bin becomes significant when it has been assigned a name and a taxon ID.  This method
     * gets the list of significant bins.
     *
     * @return the set of significant bins in this group
     */
    public List<Bin> getSignificantBins() {
        List<Bin> retVal = this.binList.stream().filter(x -> x.isSignificant()).collect(Collectors.toList());
        return retVal;
    }

    /**
     * Find the bin for a contig.
     *
     * @param contigId	ID of desired contig
     *
     * @return the bin containing the contig, or NULL if none
     */
    public Bin getContigBin(String contigId) {
        return this.contigMap.get(contigId);
    }

    /**
     * Save this bin group to a file.
     *
     * @param outFile	file in which to save the bin group
     *
     * @throws IOException
     */
    public void save(File outFile) throws IOException {
        // Convert the bin group to a json string.
        JsonObject json = this.toJson();
        String jsonString = Jsoner.serialize(json);
        // Write it out in pretty format.
        try (PrintWriter writer = new PrintWriter(outFile)) {
            writer.write(Jsoner.prettyPrint(jsonString));
        }
        log.info("Bin group saved to {}.", outFile);
    }

    /**
     * @return a JSON object for this bin group
     */
    protected JsonObject toJson() {
        JsonObject retVal = new JsonObject();
        // Form the bins into a JSON list.
        JsonArray binArray = new JsonArray();
        for (Bin bin : this.binList)
            binArray.add(bin.toJson());
        retVal.put(GroupKeys.BINS.getKey(), binArray);
        // Form the counts into a JSON object.
        JsonObject counts = new JsonObject();
        for (var counter : this.stats.counts())
            counts.put(counter.getKey(), counter.getCount());
        retVal.put(GroupKeys.COUNTS.getKey(), counts);
        // Store the file name.
        if (this.inputFile != null)
            retVal.put(GroupKeys.IN_FILE.getKey(), this.inputFile.getAbsolutePath());
        return retVal;
    }

    /**
     * Write all the unplaced contigs to a FASTA file.
     *
     * @param inFile	FASTA file containing the contigs
     * @param outFile	FASTA file to contain the unplaced contigs
     *
     * @throws IOException
     */
    public void writeUnplaced(File inFile, File outFile) throws IOException {
        // Open the two files.
        log.info("Transferring unplaced sequences from {} to {}.", inFile, outFile);
        try (FastaInputStream inStream = new FastaInputStream(inFile);
                FastaOutputStream outStream = new FastaOutputStream(outFile)) {
            int outCount = 0;
            int skipCount = 0;
            int placeCount = 0;
            // Read the input and copy to the output if warranted.
            for (Sequence seq : inStream) {
                Bin contigBin = this.contigMap.get(seq.getLabel());
                if (contigBin == null) {
                    // Here the sequence was filtered out prior to binning.
                    skipCount++;
                } else if (contigBin.isSignificant()) {
                    // Here the contig is placed in a bin.
                    placeCount++;
                } else {
                    // Here the contig should be output.
                    outStream.write(seq);
                    outCount++;
                }
            }
            log.info("{} contigs are placed, {} have been rejected, {} written to {}.", placeCount, skipCount, outCount, outFile);
        }
    }
    /**
     * Increment a statistic.
     *
     * @param statName		name of statistic to increment
     */
    public void count(String statName) {
        this.stats.count(statName);
    }

    /**
     * Increment a statistic by a specified amount.
     *
     * @param statName		name of statistic to increment
     * @param count			value to add to the statistic
     */
    public void count(String statName, int count) {
        this.stats.count(statName, count);
    }

    /**
     * @return a count
     *
     * @param statName		name of count to return
     */
    public int getCount(String statName) {
        return this.stats.getCount(statName);
    }

    /**
     * @return the number of bins in this group
     */
    public int size() {
        return this.binList.size();
    }

    /**
     * @return the input file for the bin sequences
     */
    public File getInputFile() {
        return this.inputFile;
    }

    /**
     * Specify the input file for the bin sequences.  This is protected, since it is only used for testing.
     *
     * @param inputFile 	the input file to set
     */
    protected void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }

    @Override
    public Iterator<Bin> iterator() {
        return this.binList.iterator();
    }


}
