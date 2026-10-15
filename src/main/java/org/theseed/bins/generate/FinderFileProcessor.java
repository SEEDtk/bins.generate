package org.theseed.bins.generate;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.io.LineReader;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This command creates a FinderKmer load file from a protein finder. The finder must have been previously sorted using
 * SortFinderProcessor. We first get the list of roles to form the output header. Then we do a merge load on all the
 * FASTA files and build a map of all the feature IDs we need to build the first batch of output records. These
 * IDs are then used to extract the corresponding protein sequences from the database.
 * 
 * The positional parameter is the name of the protein finder directory. The load file will be produced on the standard
 * output. It is a tab-delimited file, with headers, each line containing a genome ID followed by the protein strings 
 * for the protein finder's roles.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -o   name for the output load file (if not STDOUT)
 * -b   number of genomes per load batch (default 200)
 * 
 */
public class FinderFileProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */ 
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderFileProcessor.class);
    /** list of role IDs for the output header */
    private List<String> roleIds;
    /** protein finder manager object */
    private ProteinFinder finder;
    /** map of genome IDs to feature ID lists for the current batch */
    private Map<String, List<Set<String>>> genomeFeatures;
    /** array of IDs of the next feature ID to be read from each FASTA file */
    private String[] nextFeatures;
    /** array of line readers for each FASTA file */
    private LineReader[] readers;
    /** this will be the size for the batch feature set */
    private int batchFeatureSetSize;
    /** count of warning messages */
    private int warningCount;
    /** connection to the BV-BRC database */
    private P3CursorConnection p3;
    /** ID of the key role */
    private static final String KEY_ROLE_ID = "PhenTrnaSyntAlph";
    /** feature ID to use for end-of-file */
    private static final String EOF_FEATURE_ID = "fig|Z.peg.1";

    // COMMAND-LINE OPTIONS

    /** number of genomes per load batch */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "300", usage = "number of genomes per load batch")
    private int batchSize;

    /** name of protein-finder directory */
    @Argument(index = 0, metaVar = "FinderDir", usage = "name of protein-finder directory", required = true)
    private File finderDir;

    @Override
    protected void setReporterDefaults() {
        this.batchSize = 200;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Validate the batch size.
        if (this.batchSize <= 0)
            throw new IOException("Batch size must be positive.");
        // Create the genome batch map.
        this.genomeFeatures = new HashMap<>(this.batchSize * 2);
        // Insure the protein finder directory exists.
        if (! this.finderDir.isDirectory())
            throw new IOException("Protein finder directory " + this.finderDir + " is not found or not a directory.");
        // Load the protein finder.
        this.finder = new ProteinFinder(this.finderDir);
        // Get the role IDs.
        File roleFile = this.finder.getRoleFile();
        RoleMap finderRoles = RoleMap.load(roleFile);
        this.roleIds = new ArrayList<>(finderRoles.keySet());
        log.info("Loaded {} role IDs from {}", this.roleIds.size(), roleFile);
        // The key role ID must be first in the role list. Every genome has an instance of the
        // key role, so we can use it to simplify the merge process.
        if (! this.roleIds.contains(KEY_ROLE_ID))
            throw new IOException("Key role ID " + KEY_ROLE_ID + " not found in role file.");
        this.roleIds.remove(KEY_ROLE_ID);
            this.roleIds.add(0, KEY_ROLE_ID);
        // Connect to the BV-BRC database.
        this.p3 = new P3CursorConnection();
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write out the header line.
        writer.print("genome_id");
        for (String roleId : this.roleIds) {
            writer.print("\t");
            writer.print(roleId);
        }
        writer.println();
        // This will count the number of missing proteins in the database.
        this.warningCount = 0;
        // This will count the number of genomes output.
        int genomeCount = 0;
        // This will count the number of missing proteins.
        int missingCount = 0;
        // This will count the number of proteins read.
        int proteinCount = 0;
        // This will count the number of batches processed.
        int batchCount = 0;
        // Now we create an array of line readers for the FASTA files.
        this.readers = new LineReader[this.roleIds.size()];
        // Compute a good initial size for the batch feature set. A full set of roles for each
        // genome in the batch.
        this.batchFeatureSetSize = this.roleIds.size() * this.batchSize * 3 / 2 + 1;
        // FASTA files require a read-ahead logic, so we create an array of
        // buffers for the next feature ID from each FASTA file.
        this.nextFeatures = new String[this.readers.length];
        try {
            this.openReaders();
            // Process each genome ID in the key role's file (which is position 0 in the array).
            while (! this.nextFeatures[0].equals(EOF_FEATURE_ID)) {
                // Insure there is room for this next genome.
                if (this.genomeFeatures.size() >= this.batchSize) {
                    batchCount++;
                    log.info("Processing batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
                    this.processBatch(writer);
                    this.genomeFeatures.clear();
                }
                // Get the current genome's ID from the feature ID.
                String genomeId = Feature.genomeOf(this.nextFeatures[0]);
                genomeCount++;
                // Put all of its features into a list.
                List<Set<String>> features = new ArrayList<>(readers.length);
                int emptyCount = 0;
                for (int i = 0; i < this.readers.length; i++) {
                    // Initialize the feature slot for this reader. So long as we have a
                    // matching genome ID, we keep reading features from this reader. It is
                    // rare to have multiple feature IDs for a universal role, but sadly it
                    // does happen.
                    Set<String> featureSet = new TreeSet<>();
                    while (Feature.genomeOf(this.nextFeatures[i]).equals(genomeId)) {
                        featureSet.add(this.nextFeatures[i]);
                        proteinCount++;
                        this.readNext(i);
                    }
                    if (featureSet.isEmpty()) {
                        missingCount++;
                        emptyCount++;
                    }
                    features.add(featureSet);
                }
                if (emptyCount > 3)
                    log.warn("Genome {} is missing {} features.", genomeId, emptyCount);
                // Add the genome to the batch.
                this.genomeFeatures.put(genomeId, features);
            }
            // Process any remaining genomes in the batch.
            if (! this.genomeFeatures.isEmpty()) {
                batchCount++;
                log.info("Processing final batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
                this.processBatch(writer);
            }
            log.info("All done. {} genomes processed. {} proteins read. {} missing proteins.", genomeCount, proteinCount, missingCount);
        } finally {
            int closeCount = 0;
            for (LineReader reader : readers) {
                if (reader != null) {
                    reader.close();
                    closeCount++;
                }
            }
            log.info("Closed {} FASTA readers.", closeCount);
        }
        if (this.warningCount > 0)
            log.warn("There were {} warnings during processing.", this.warningCount);
        log.info("Processed {} genomes with {} missing proteins and {} proteins read.", genomeCount, missingCount, proteinCount);
    }

    /**
     * Open the FASTA files corresponding to the role IDs and store the readers in the array.
     * 
     * @throws IOException if a FASTA file cannot be opened
     */
    private void openReaders() throws IOException {
        Map<String, File> fastaFiles = this.finder.getFastas();
        for (int i = 0; i < this.readers.length; i++) {
            String roleId = this.roleIds.get(i);
            File fastaFile = fastaFiles.get(roleId);
            log.info("Role {} feature IDs will be read from {}.", roleId, fastaFile);
            this.readers[i] = new LineReader(fastaFile);
            // The first line will contain the first feature ID, so we read it now. If a FASTA file is empty,
            // this will throw an exception. This is not a problem since a universal role that never appears
            // is not allowed.
            String line = this.readers[i].next();
            if (! line.startsWith(">"))
                throw new IOException("Expected FASTA header line for role " + roleId + " in file " + fastaFile);
            // The feature ID is everything after the '>' up to the first space.
            this.nextFeatures[i] = StringUtils.substringBefore(line.substring(1), " ");
        }
    }

    /**
     * Read the next feature ID from a FASTA file. Return the end-of-file feature ID if the file is exhausted.
     *
     * @param i	index of the FASTA reader
     * 
     * @return the next feature ID from the FASTA file, or the end-of-file feature ID if the file is exhausted
     * 
     * @throws IOException if an I/O error occurs
     */
    private void readNext(int i) throws IOException {
        boolean found = false;
        while (this.readers[i].hasNext() && ! found) {
            String line = this.readers[i].next();
            if (line.startsWith(">")) {
                this.nextFeatures[i] = StringUtils.substringBefore(line.substring(1), " ");
                found = true;
            }
        }
        if (! found)
            this.nextFeatures[i] = EOF_FEATURE_ID;
    }

    /**
     * Process the current batch of genomes.
     * 
     * @param writer the writer for writing the batch output lines
     * 
     * @throws IOException if a database error occurs
     */
    private void processBatch(PrintWriter writer) throws IOException {
        // Gather all the features for the current batch of genomes into a set.
        Set<String> batchFeatures = new HashSet<>(this.batchFeatureSetSize);
        this.genomeFeatures.values().stream().flatMap(x -> x.stream()).forEach(x -> batchFeatures.addAll(x));
        log.info("Gathered {} features for the current batch.", batchFeatures.size());
        Map<String,String> featureSeqs = this.getFeatureSequences(batchFeatures);
        // Now we loop through the genomes and write the proteins.
        for (var genomeEntry : this.genomeFeatures.entrySet()) {
            String genomeId = genomeEntry.getKey();
            writer.print(genomeId);
            // The tricky part here is that sometimes there is more than one feature for a genome/role pair,
            // in which case we output the proteins comma-delimited.
            List<Set<String>> featureSets = genomeEntry.getValue();
            for (Set<String> features : featureSets) {
                writer.print("\t");
                if (! features.isEmpty()) {
                    List<String> proteins = new ArrayList<>(2);
                    for (String feature : features) {
                        String protein = featureSeqs.get(feature);
                        if (protein == null) {
                            log.warn("Protein sequence not found for feature {} in genome {}", feature, genomeId);
                            this.warningCount++;
                        } else
                            proteins.add(protein);
                    }
                    if (! proteins.isEmpty())
                        writer.print(String.join(",", proteins));
                }
            }
            // Terminate the current genome line.
            writer.println();
        }
        // Insure the batch is set to the output.
        writer.flush();
    }

    /**
     * Get all the feature sequences from the batch of features.
     * 
     * @param batchFeatures     set of feature IDs for which sequences are to be retrieved
     * 
     * @return a map of feature IDs to their corresponding sequences
     * 
     * @throws IOException if a database error occurs
     */
    private Map<String, String> getFeatureSequences(Set<String> batchFeatures) throws IOException {
        Map<String, String> retVal = new HashMap<>(this.batchFeatureSetSize);
        this.p3.getRecords("feature", P3CursorConnection.MAX_LIMIT, batchFeatures.size(), "patric_id", batchFeatures, "patric_id,aa_sequence",
            Collections.emptySet(), x -> retVal.put(KeyBuffer.getString(x, "patric_id"), KeyBuffer.getString(x, "aa_sequence")));
        return retVal;
    }

}
