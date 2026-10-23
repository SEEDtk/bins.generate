package org.theseed.bins.generate;

import java.io.File;
import java.io.FileWriter;
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

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseInputProcessor;

/**
 * This command creates a FinderKmer load file from the evaluation master database. This is a directory containing
 * compressed files for all the genomes in the BV-BRC. The FinderKmer load file contains the ID of each good genome
 * along with the protein sequences of all its key universal proteins. We get the IDs of the universal proteins
 * from the master database and then extract the actual protein sequences from the main BV-BRC database. The output
 * file is tab-delimited with headers.
 * 
 * The positional parameters are the name of the evaluation master database, the definition file for the roles to use,
 * and the name of the output file.
 * 
 * The standard input should be the list of good genomes from the BV-BRC; however, all that is really needed is a 
 * tab-delimited file with genome IDs in the first column.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -i   name of the input file containing the list of good genomes (if not STDIN)
 * -b   number of genomes per load batch (default 200)
 * -c   index (1-based) or name of the input column containing the genome IDs (default "1")
 * 
 * --resume     if specified, the output file will be read and processing will continue from where it left off
 * 
 */
public class KmerLoadFileProcessor extends BaseInputProcessor {

    // FIELDS
    /** logging facility */ 
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(KmerLoadFileProcessor.class);
    /** list of role IDs for the output header */
    private List<String> roleIds;
    /** this will be the size for the batch feature set */
    private int batchFeatureSetSize;
    /** count of warning messages */
    private int warningCount;
    /** connection to the BV-BRC database */
    private P3CursorConnection p3;
    /** map of genome IDs to feature IDs (in order) for the current batch */
    private Map<String, List<Set<String>>> genomeFeatures;
    /** list of genome IDs in input order */
    private List<String> genomeIds;
    /** map of role IDs to array positions */
    private Map<String, Integer> rolePositions;
    /** genome source for the evaluation master database */
    private GenomeSource genomes;
    /** input column index for genome IDs */
    private int genomeIdColIdx;
    /** map of role IDs to Role definitions */
    private RoleMap roleMap;
    /** ID of the last processed genome (for resuming) */
    private String lastProcessedId;
    /** number of genomes skipped due to resume */
    private int skipped;

    // COMMAND-LINE OPTIONS

    /** number of genomes per load batch */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "300", usage = "number of genomes per load batch")
    private int batchSize;

    /** index (1-based) pr name of the input column containing the genome IDs */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "1", usage = "index (1-based) or name of the input column containing the genome IDs")
    private String genomeIdCol;

    /** if specified, we are resuming from a previous run */
    @Option(name = "--resume", usage = "if specified, we are resuming from a previous run")
    private boolean resume;

    /** name of evaluation master database */
    @Argument(index = 0, metaVar = "P3Master", usage = "name of evaluation master database", required = true)
    private File evalMasterDir;
    
    /** definition file for the roles to use */
    @Argument(index = 1, metaVar = "roles.for.finder", usage = "name of the definition file for the roles to use", required = true)
    private File roleFile;

    /** output file for the report */
    @Argument(index = 2, metaVar = "loadFile.tbl", usage = "output file for the report", required = true)
    private File outFile;

    @Override
    protected void setReaderDefaults() {
        this.batchSize = 200;
        this.genomeIdCol = "1";
        this.resume = false;
    }

    @Override
    protected void validateReaderParms() throws IOException, ParseFailureException {
        // Validate the batch size.
        if (this.batchSize <= 0)
            throw new IOException("Batch size must be positive.");
        // Create the genome batch map and list.
        this.genomeFeatures = new HashMap<>(this.batchSize * 2);
        this.genomeIds = new ArrayList<>(this.batchSize);
        // Insure the protein finder directory exists.
        if (! this.evalMasterDir.isDirectory())
            throw new IOException("Evaluation master database directory " + this.evalMasterDir + " is not found or not a directory.");
        // Load the protein finder.
        this.genomes = GenomeSource.Type.MASTER.create(this.evalMasterDir);
        // Get the role IDs.
        this.roleMap = RoleMap.load(this.roleFile);
        this.roleIds = new ArrayList<>(roleMap.keySet());
        this.rolePositions = new HashMap<>();
        for (int i = 0; i < this.roleIds.size(); i++)
            this.rolePositions.put(this.roleIds.get(i), i);
        log.info("Loaded {} role IDs from {}", this.roleIds.size(), roleFile);
        // Check if we are resuming from a previous run.
        this.lastProcessedId = null;
        if (this.resume) {
            log.info("Resuming from a previous run.");
            if (! this.outFile.canRead())
                throw new IOException("Cannot resume because output file " + this.outFile + " is not readable.");
            // Find the last processed line in the output file.
            try (TabbedLineReader reader = new TabbedLineReader(this.outFile)) {
                this.lastProcessedId = null;
                this.skipped = 0;
                for (var line : reader) {
                    this.lastProcessedId = line.get(this.genomeIdColIdx);
                    this.skipped++;
                }
                log.info("{} genomes were previously output. Resuming after genome ID {}.", this.skipped, this.lastProcessedId);
            }
        }
        // Connect to the BV-BRC database.
        this.p3 = new P3CursorConnection();
    }

    @Override
    protected void validateReaderInput(TabbedLineReader input) throws IOException {
        // This will throw an IOException if the genome ID column is not found.
        this.genomeIdColIdx = input.findField(this.genomeIdCol);
    }

    @Override 
    protected void runReader(TabbedLineReader input) throws Exception {
        // Open the output writer.
        try (PrintWriter writer = this.openOutputWriter()) {
            // This will count the number of missing proteins in the database.
            this.warningCount = 0;
            // This will count the number of missing genomes in the master database.
            int missingGenomes = 0;
            // This will count the number of genomes output.
            int genomeCount = 0;
            // This will count the number of missing proteins.
            int missingCount = 0;
            // This will count the number of proteins read.
            int proteinCount = 0;
            // This will count the number of batches processed.
            int batchCount = 0;
            // This will count the number of genomes skipped due to resuming.
            int skippedGenomes = 0;
            // Loop through the input, loading each genome from the master database.
            for (var line : input) {
                // Get the new genome ID.
                String genomeId = line.get(this.genomeIdColIdx);
                // Check for a resume. If we have a last-processed genome, skip until we get past it.
                if (this.lastProcessedId != null) {
                    if (! genomeId.equals(this.lastProcessedId))
                        skippedGenomes++;
                    else {
                        this.lastProcessedId = null;
                        log.info("Resuming processing after skipping {} genomes.", skippedGenomes);
                    }
                } else {
                    // Here we are ready to process. Insure there is room in the batch for a new genome.
                    if (this.genomeFeatures.size() >= this.batchSize) {
                        batchCount++;
                        log.info("Processing batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
                        this.processBatch(writer);
                        this.genomeFeatures.clear();
                        this.genomeIds.clear();
                        log.info("{} total genomes in output file.", genomeCount + this.skipped);
                    }
                    // Denote we're processing a new genome.
                    genomeCount++;
                    Genome genome = this.genomes.getGenome(genomeId, P3Genome.Details.STRUCTURE_ONLY);
                    if (genome == null) {
                        log.warn("Genome {} not found in the database.", genomeId);
                        missingGenomes++;
                    } else {
                        // This genome exists, so add it to the batch.
                        this.genomeIds.add(genomeId);
                        // Extract the features we want from this genome and store them in the batch.
                        for (Feature feat : genome.getFeatures()) {
                            List<Role> roles = feat.getUsefulRoles(this.roleMap);
                            for (Role role : roles) {
                                String roleId = role.getId();
                                if (roleId != null && this.rolePositions.containsKey(roleId)) {
                                    int rolePos = this.rolePositions.get(roleId);
                                    List<Set<String>> genomeFeatureIds = this.genomeFeatures.computeIfAbsent(genomeId, k -> this.newSourList());
                                    Set<String> roleFeatures = genomeFeatureIds.get(rolePos);
                                    roleFeatures.add(feat.getId());
                                    proteinCount++;
                                }
                            }
                        }
                        int emptyCount = 0;
                        for (Set<String> roleFeatures : this.genomeFeatures.get(genomeId)) {
                            if (roleFeatures.isEmpty())
                                emptyCount++;
                        }
                        if (emptyCount > 3)
                            log.warn("Genome {} is missing {} features.", genomeId, emptyCount);
                        missingCount += emptyCount;
                    }
                }
            }
            // Process any remaining genomes in the batch.
            if (! this.genomeFeatures.isEmpty()) {
                batchCount++;
                log.info("Processing final batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
                this.processBatch(writer);
                log.info("{} total genomes in output file.", genomeCount + this.skipped);
            }
            log.info("All done. {} genomes processed. {} proteins read. {} missing proteins.", genomeCount, proteinCount, missingCount);
            if (this.warningCount > 0)
                log.warn("There were {} missing proteins in the database.", this.warningCount);
            if (missingGenomes > 0)
                log.warn("There were {} genomes missing from the database.", missingGenomes);
            log.info("Processed {} genomes with {} missing proteins and {} proteins read.", genomeCount, missingCount, proteinCount);
        }
    }

    /**
     * Open the output file for writing. We open for append if we are resuming. For a fresh run,
     * we open for writing and write a header line.
     * 
     * @return a PrintWriter for writing the output file
     * 
     * @throws IOException if the output file cannot be opened for writing
     */
    private PrintWriter openOutputWriter() throws IOException {
        PrintWriter retVal;
        if (this.resume) {
            log.info("Appending new genomes to existing output file.");
            retVal = new PrintWriter(new FileWriter(this.outFile, true));
        } else {
            retVal = new PrintWriter(new FileWriter(this.outFile));
            // Write the header line for a fresh run.
            retVal.print("genome_id");
            this.roleIds.stream().forEach(roleId -> retVal.print("\t" + roleId));
        }
        return retVal;
    }

    /**
     * Create an empty skeleton for a genome's feature list.
     * 
     * @return a list of empty sets, one for each role position
     */
    private List<Set<String>> newSourList() {
        List<Set<String>> retVal = new ArrayList<>();
        for (int i = 0; i < this.rolePositions.size(); i++)
            retVal.add(new TreeSet<>());
        return retVal;
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
        for (String genomeId : this.genomeIds) {
            writer.print(genomeId);
            // The tricky part here is that sometimes there is more than one feature for a genome/role pair,
            // in which case we output the proteins comma-delimited.
            List<Set<String>> featureSets = this.genomeFeatures.get(genomeId);
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
