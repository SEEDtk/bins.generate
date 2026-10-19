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
import org.theseed.utils.BasePipeProcessor;

/**
 * This command creates a FinderKmer load file from the evaluation master database. This is a directory containing
 * compressed files for all the genomes in the BV-BRC. The FinderKmer load file contains the ID of each good genome
 * along with the protein sequences of all its key universal proteins. We get the IDs of the universal proteins
 * from the master database and then extract the actual protein sequences from the main BV-BRC database. The output
 * file is tab-delimited with headers.
 * 
 * The positional parameters are the name of the evaluation master database and the definition file for the roles to use. 
 * The standard input should be the list of good genomes from the BV-BRC; however, all that is really needed is a 
 * tab-delimited file with genome IDs in the first column. The output load file will be written to the standard output.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -i   name of the input file containing the list of good genomes (if not STDIN)
 * -o   name for the output load file (if not STDOUT)
 * -b   number of genomes per load batch (default 200)
 * -c   index (1-based) or name of the input column containing the genome IDs (default "1")
 * 
 */
public class KmerLoadFileProcessor extends BasePipeProcessor {

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
    /** map of role IDs to array positions */
    private Map<String, Integer> rolePositions;
    /** genome source for the evaluation master database */
    private GenomeSource genomes;
    /** input column index for genome IDs */
    private int genomeIdColIdx;
    /** map of role IDs to Role definitions */
    private RoleMap roleMap;

    // COMMAND-LINE OPTIONS

    /** number of genomes per load batch */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "300", usage = "number of genomes per load batch")
    private int batchSize;

    /** index (1-based) pr name of the input column containing the genome IDs */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "1", usage = "index (1-based) or name of the input column containing the genome IDs")
    private String genomeIdCol;

    /** name of evaluation master database */
    @Argument(index = 0, metaVar = "P3Master", usage = "name of evaluation master database", required = true)
    private File evalMasterDir;
    
    /** definition file for the roles to use */
    @Argument(index = 1, metaVar = "roles.for.finder", usage = "name of the definition file for the roles to use", required = true)
    private File roleFile;

    @Override
    protected void setPipeDefaults() {
        this.batchSize = 200;
        this.genomeIdCol = "1";
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Validate the batch size.
        if (this.batchSize <= 0)
            throw new IOException("Batch size must be positive.");
        // Create the genome batch map.
        this.genomeFeatures = new HashMap<>(this.batchSize * 2);
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
        // Connect to the BV-BRC database.
        this.p3 = new P3CursorConnection();
    }

    @Override
    protected void validatePipeInput(TabbedLineReader input) throws IOException {
        // This will throw an IOException if the genome ID column is not found.
        this.genomeIdColIdx = input.findField(this.genomeIdCol);
    }

    @Override 
    protected void runPipeline(TabbedLineReader input, PrintWriter writer) throws Exception {
        // Write out the header line.
        writer.print("genome_id");
        for (String roleId : this.roleIds) {
            writer.print("\t");
            writer.print(roleId);
        }
        writer.println();
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
        // Loop through the input, loading each genome from the master database.
        for (var line : input) {
            // Insure there is room in the batch for a new genome.
            if (this.genomeFeatures.size() >= this.batchSize) {
                batchCount++;
                log.info("Processing batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
                this.processBatch(writer);
                this.genomeFeatures.clear();
            }
            // Get the new genome ID.
            String genomeId = line.get(this.genomeIdColIdx);
            genomeCount++;
            Genome genome = this.genomes.getGenome(genomeId, P3Genome.Details.STRUCTURE_ONLY);
            if (genome == null) {
                log.warn("Genome {} not found in the database.", genomeId);
                missingGenomes++;
            } else {
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
            }
        }
        // Process any remaining genomes in the batch.
        if (! this.genomeFeatures.isEmpty()) {
            batchCount++;
            log.info("Processing final batch {} of {} genomes.", batchCount, this.genomeFeatures.size());
            this.processBatch(writer);
        }
            log.info("All done. {} genomes processed. {} proteins read. {} missing proteins.", genomeCount, proteinCount, missingCount);
        if (this.warningCount > 0)
            log.warn("There were {} missing proteins in the database.", this.warningCount);
        if (missingGenomes > 0)
            log.warn("There were {} genomes missing from the database.", missingGenomes);
        log.info("Processed {} genomes with {} missing proteins and {} proteins read.", genomeCount, missingCount, proteinCount);
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
