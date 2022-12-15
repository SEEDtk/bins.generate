/**
 *
 */
package org.theseed.bins;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.NavigableSet;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object describes a bin.  A bin has a computed coverage and various bits of information
 * that can be used to tell RAST how to annotate it.  It is also associated with a list of contigs.  A bin can
 * be serialized to and from JSON format.
 *
 * @author Bruce Parrello
 *
 */
public class Bin implements Comparable<Bin> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(Bin.class);
    /** name of the bin, usually derived from the recommended taxon */
    private String name;
    /** taxon ID of the bin */
    private int taxonID;
    /** domain of the bin (Bacteria, Archaea) */
    private String domain;
    /** genetic code of the bin */
    private int gc;
    /** reference genome ID list, ordered from best to worst */
    private List<String> refGenomes;
    /** set of contigs */
    private NavigableSet<Member> contigs;
    /** status of the bin (not saved to JSON) */
    private Status status;
    /** output file for this bin */
    private File outFile;
    /** output stream for this bin (not saved to JSON) */
    private FastaOutputStream outStream;
    /** constant for an empty bin */
    private static final Set<JsonArray> noContigs = Collections.emptySet();
    /** constant for no reference genomes */
    private static final List<String> noRefGenomes = Collections.emptyList();
    /** default taxon ID */
    private static final int DEFAULT_TAX_ID = 0;
    /** default domain */
    private static final String DEFAULT_DOMAIN = "Bacteria";
    /** default genetic code */
    private static final int DEFAULT_GC = 11;

    /**
     * This enum indicates what should be done with a contig bin during import.
     */
    public static enum Status {
        /** discard the contig */
        BAD,
        /** use the contig for binning only */
        NORMAL,
        /** use the contig for binning and for the SOUR protein search */
        SEED_USABLE;

        /**
         * @return a status name for use in statistical counters
         */
        public String toString() {
            return this.name().toLowerCase() + "_contig";
        }
    }


    /**
     * This enum is used to convert the bin to and from JSON.
     */
    protected static enum BinKeys implements JsonKey {
        NAME("unknown species bin"),
        GC(DEFAULT_GC),
        DOMAIN(DEFAULT_DOMAIN),
        TAXON_ID(DEFAULT_TAX_ID),
        REF_GENOME(noRefGenomes),
        FILE_NAME(null),
        CONTIGS(noContigs);

        /** default value for key */
        private Object mValue;

        private BinKeys(Object value) {
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
     * This class represents a contig member of a bin.  In JSON form, it is a 3-element
     * list:  contig ID, length, coverage.  Members are ordered from longest to shortest,
     * then best coverage to worst, then by contig ID.
     */
    public static class Member implements Comparable<Member> {

        /** ID of the contig */
        private String contigId;
        /** length of the contig */
        private int len;
        /** coverage of the contig */
        private double coverage;

        /**
         * Create a new contig member.
         *
         * @param contigId		ID of the contig
         * @param len			length of the contig
         * @param coverage		coverage of the contig
         */
        public Member(String contigId, int len, double coverage) {
            this.contigId = contigId;
            this.len = len;
            this.coverage = coverage;
        }

        /**
         * Create a contig member from a JSON array.
         *
         * @param json	JSON array from which to build the member
         */
        public Member(JsonArray json) {
            this.contigId = json.getString(0);
            this.len = json.getInteger(1);
            this.coverage = json.getDouble(2);
        }

        /**
         * @return a JSON array representing this member
         */
        public JsonArray toJson() {
            JsonArray retVal = new JsonArray();
            retVal.addChain(this.contigId).addChain(this.len).addChain(this.coverage);
            return retVal;
        }

        /**
         * @return the ID of this contig
         */
        public String getContigId() {
            return this.contigId;
        }

        /**
         * @return the length of this contig
         */
        public int getLen() {
            return this.len;
        }

        /**
         * @return the coverage level of this contig
         */
        public double getCoverage() {
            return this.coverage;
        }

        @Override
        public int compareTo(Member o) {
            int retVal = o.len - this.len;
            if (retVal == 0) {
                retVal = Double.compare(o.coverage, this.coverage);
                if (retVal == 0)
                    retVal = this.contigId.compareTo(o.contigId);
            }
            return retVal;
        }


        @Override
        public int hashCode() {
            final int prime = 31;
            int result = prime + ((this.contigId == null) ? 0 : this.contigId.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Member)) {
                return false;
            }
            Member other = (Member) obj;
            if (this.contigId == null) {
                if (other.contigId != null) {
                    return false;
                }
            } else if (!this.contigId.equals(other.contigId)) {
                return false;
            }
            return true;
        }

    }

    /**
     * Create a bin from a single contig.
     *
     * @param contigId		ID of the contig
     * @param len			length of the contig
     * @param covg			coverage of the contig
     */
    public Bin(String contigId, int len, double covg) {
        this.contigs = new TreeSet<Member>();
        this.contigs.add(new Member(contigId, len, covg));
        // The initial name is the contig ID.
        this.name = contigId;
        // Default the other stuff.
        this.taxonID = DEFAULT_TAX_ID;
        this.gc = DEFAULT_GC;
        this.domain = DEFAULT_DOMAIN;
        this.refGenomes = new ArrayList<String>(5);
        this.outFile = null;
        this.outStream = null;
    }

    /**
     * Construct this bin from a JSON object.
     *
     * @param json	JSON object describing the bin
     */
    public Bin(JsonObject json) {
        this.taxonID = json.getIntegerOrDefault(BinKeys.TAXON_ID);
        this.refGenomes = json.getCollectionOrDefault(BinKeys.REF_GENOME);
        this.gc = json.getIntegerOrDefault(BinKeys.GC);
        this.name = json.getStringOrDefault(BinKeys.NAME);
        this.domain = json.getStringOrDefault(BinKeys.DOMAIN);
        // Get the file name (if any).
        String fileString = json.getStringOrDefault(BinKeys.FILE_NAME);
        if (fileString == null)
            this.outFile = null;
        else
            this.outFile = new File(fileString);
        // Read in the contigs.
        Collection<JsonArray> contigList = json.getCollectionOrDefault(BinKeys.CONTIGS);
        this.contigs = new TreeSet<Member>();
        for (JsonArray contig : contigList) {
            this.contigs.add(new Member(contig));
        }
    }

    /**
     * @return a JSON object describing this bin
     */
    public JsonObject toJson() {
        JsonObject retVal = new JsonObject();
        retVal.put(BinKeys.DOMAIN.getKey(), this.domain);
        retVal.put(BinKeys.GC.getKey(), this.gc);
        retVal.put(BinKeys.NAME.getKey(), this.name);
        retVal.put(BinKeys.TAXON_ID.getKey(), this.taxonID);
        if (this.outFile != null)
            retVal.put(BinKeys.FILE_NAME.getKey(), this.outFile.getAbsolutePath());
        // Form the reference genomes into a JSON array.
        JsonArray refGenomeList = new JsonArray(this.refGenomes);
        retVal.put(BinKeys.REF_GENOME.getKey(), refGenomeList);
        // Form the members into a JSON array.
        JsonArray members = new JsonArray();
        for (Member member : this.contigs)
            members.add(member.toJson());
        retVal.put(BinKeys.CONTIGS.getKey(), members);
        return retVal;
    }

    /**
     * @return the domain for this bin
     */
    public String getDomain() {
        return this.domain;
    }

    /**
     * @return the name of this bin
     */
    public String getName() {
        return this.name;
    }

    /**
     * @return the recommended taxonomic grouping for this bin
     */
    public int getTaxonID() {
        return this.taxonID;
    }

    /**
     * @return the genetic code of this bin
     */
    public int getGc() {
        return this.gc;
    }

    /**
     * @return the contigs of this bin
     */
    public List<String> getContigs() {
        return this.contigs.stream().map(x -> x.getContigId()).collect(Collectors.toList());
    }

    /**
     * Specify the output file for this bin.
     *
     * @param outFile	name of the output file
     *
     * @return the open output stream for the bin
     *
     * @throws IOException
     */
    public FastaOutputStream setOutFile(File outFileName) throws IOException {
        log.info("Bin {} will be written to {}.", this.name, outFileName);
        this.outFile = outFileName;
        this.outStream = new FastaOutputStream(outFileName);
        return this.outStream;
    }

    /**
     * @return TRUE if this bin has an open output stream, else FALSE
     */
    protected boolean isOpen() {
        return this.outStream != null;
    }

    /**
     * Write a FASTA sequence to this bin's output stream.
     *
     * @param seq	sequence to write
     *
     * @throws IOException
     */
    public void writeSequence(Sequence seq) throws IOException {
        this.outStream.write(seq);
    }

    /**
     * Close the output file for this bin.
     */
    public void close() {
        if (this.outStream != null) {
            this.outStream.close();
            this.outStream = null;
        }
    }

    /**
     * @return the name of the output file to which this bin was written, or NULL if it has never been written
     */
    public File getOutFile() {
        return this.outFile;
    }

    /**
     * @return the coverage of this bin
     */
    public double getCoverage() {
        double retVal = 0.0;
        int totLen = 0;
        for (var member : this.contigs) {
            int len = member.getLen();
            retVal += member.getCoverage() * len;
            totLen += len;
        }
        if (retVal > 0.0)
            retVal /= totLen;
        return retVal;
    }

    /**
     * @return the total length of this bin
     */
    public int getLen() {
        int retVal = this.contigs.stream().mapToInt(x -> x.getLen()).sum();
        return retVal;
    }

    /**
     * Incorporate a species recommendation into this bin.
     *
     * @param recommendation	recommendation from the protein finder
     * @param suffix			suffix to add to the bin name
     * @param refGenome			reference genome object
     */
    public void setTaxInfo(ProteinFinder.DnaHit recommendation, String suffix, Genome refGenome) {
        this.taxonID = recommendation.getSpeciesId();
        this.name = recommendation.getName();
        this.refGenomes.clear();
        this.refGenomes.add(refGenome.getId());
        if (! StringUtils.isBlank(suffix))
            this.name += " " + suffix;
        this.gc = refGenome.getGeneticCode();
        this.domain = refGenome.getDomain();
    }

    /**
     * Add a reference genome to this bin.
     *
     * @param refGenomeId	ID of the reference genome to add
     */
    public void addRefGenome(String refGenomeId) {
        // We can't use a set for the reference-genome list because the order matters,
        // so we have to check to prevent duplicates.
        if (! this.refGenomes.contains(refGenomeId))
            this.refGenomes.add(refGenomeId);
    }

    /**
     * @return TRUE if this bin has been assigned a taxon ID for output
     */
    public boolean isSignificant() {
        return this.taxonID > 0;
    }

    /**
     * @return the ID of the primary reference genome (or NULL if there is none)
     */
    public String getRefGenome() {
        String retVal;
        if (this.refGenomes.isEmpty())
            retVal = null;
        else
            retVal = this.refGenomes.get(0);
        return retVal;
    }

    /**
     * @return the IDs of all reference genomes in a list, with the primary first
     */
    public List<String> getAllRefGenomes() {
        return this.refGenomes;
    }

    /**
     * Merge another bin's data into this one.
     *
     * @param other		other bin to merge
     */
    protected void merge(Bin other) {
        this.contigs.addAll(other.contigs);
        this.refGenomes.addAll(other.refGenomes);
    }

    /**
     * This is a stronger relationship than "equals", since it checks all the derived fields.
     *
     * @param other		other bin to check
     *
     * @return TRUE if the other bin has the same fields as this one
     */
    public boolean isClone(Bin other) {
        boolean retVal = this.domain.equals(other.domain) && this.gc == other.gc && this.name.equals(other.name)
                && this.refGenomes.equals(other.refGenomes) && this.taxonID == other.taxonID
                && this.contigs.size() == other.contigs.size();
        if (retVal) {
            // Verify both bins have the same contigs.
            var otherContigs = other.contigs;
            var iter = this.contigs.iterator();
            while (iter.hasNext() && retVal) {
                Member contig = iter.next();
                Member otherContig = otherContigs.ceiling(contig);
                retVal = contig.contigId.equals(otherContig.contigId) && contig.len == otherContig.len
                        && contig.coverage == otherContig.coverage;
            }
        }
        if (retVal) {
            // Finally, check the file names.
            if (this.outFile == null)
                retVal = (other.outFile == null);
            else
                retVal = (this.outFile.getAbsolutePath().equals(other.outFile.getAbsolutePath()));
        }
        return retVal;
    }

    /**
     * @return the load disposition of this bin
     */
    public Status getStatus() {
        return this.status;
    }

    /**
     * Specify what to do with the bin during loading.
     *
     * @param status 	the disposition during load for this bin
     */
    public void setStatus(Status status) {
        this.status = status;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.name == null) ? 0 : this.name.hashCode());
        result = prime * result + this.taxonID;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof Bin)) {
            return false;
        }
        Bin other = (Bin) obj;
        if (this.name == null) {
            if (other.name != null) {
                return false;
            }
        } else if (!this.name.equals(other.name)) {
            return false;
        }
        if (this.taxonID != other.taxonID) {
            return false;
        }
        return true;
    }

    @Override
    public int compareTo(Bin o) {
        int retVal = this.taxonID - o.taxonID;
        if (retVal == 0) {
            retVal = this.name.compareTo(o.name);
        }
        return retVal;
    }

    /**
     * Specify that this bin is the virtual bin to contain the unplaced contigs.
     *
     * @param binFile	output file for unplaced contigs
     */
    public void setVirtual(File binFile) {
        this.name = "Residual Contigs";
        this.outFile = binFile;
    }

}
