/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This command builds a protein finder.  Building the SOUR protein file is fairly fast, but the reference-genome
 * DNA files can be very expensive.  If files already exist in the output directory, they will not be rebuilt.
 * Use "--clear" to force a rebuild.
 *
 * The positional parameters are the name of the output directory, the name of the role definition file to use, the genome
 * source (file or directory) containing the genomes to use for the SOUR proteins, the name of a tab-delimited file
 * containing the IDs of the eligible reference genomes in its first column.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --source		type of genome source for the SOUR protein genomes
 * --clear		erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class BuildFinderProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(BuildFinderProcessor.class);
    /** protein finder being created */
    private ProteinFinder finder;
    /** genome source */
    private GenomeSource genomes;
    /** map of genome IDs to species IDs for the eligible reference genomes */
    private Map<String, Integer> refMap;

    // COMMAND-LINE OPTIONS

    /** genome source type */
    @Option(name = "--source", usage = "type of genome source for SOUR proteins")
    private GenomeSource.Type sourceType;

    /** if specified, the output directory will be erased before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** name of the output directory */
    @Argument(index = 0, metaVar = "outDir", usage = "name of the output directory in which to create the finder")
    private File outDir;

    /** name of the role definition file */
    @Argument(index = 1, metaVar = "roles.for.proteins", usage = "name of the role definition file for the SOUR proteins")
    private File roleFile;

    /** name of the genome source for the SOUR protein genomes */
    @Argument(index = 2, metaVar = "genomeDir", usage = "name of the genome source containing the genomes to use for the SOUR protein")
    private File genomeDir;

    /** name of the input file for the reference genome IDs */
    @Argument(index = 3, metaVar = "patric.good.tbl", usage = "name of the file containing the eligible reference-genome IDs")
    private File refListFile;

    @Override
    protected void setDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.clearFlag = false;
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        // Validate the reference genome input file.
        if (! this.refListFile.canRead())
            throw new FileNotFoundException("Reference-genome list file " + this.refListFile + " is not found or unreadable.");
        this.refMap = new HashMap<>(100000);
        try (TabbedLineReader refStream = new TabbedLineReader(this.refListFile)) {
            int genomeColIdx = refStream.findField("genome_id");
            int speciesColIdx = refStream.findField("species");
            for (var line : refStream) {
                String genomeId = line.get(genomeColIdx);
                int speciesId = line.getInt(speciesColIdx);
                if (speciesId > 0)
                    this.refMap.put(genomeId, speciesId);
            }
        }
        log.info("{} reference-genome IDs read from {}.", this.refMap.size(), this.refListFile);
        // Insure the role file is readable.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role file " + this.roleFile + " is not found or unreadable.");
        // Connect to the SOUR protein genome source.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("SOUR protein genome source " + this.genomeDir + " is not found.");
        this.genomes = this.sourceType.create(this.genomeDir);
        log.info("{} SOUR protein genomes found in {}.", this.genomes.size(), this.genomeDir);
        // Clear the output directory if requested.
        if (this.outDir.isDirectory() && this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Protein finder will be created in {}.", this.outDir);
        // Initialize the protein finder.
        this.finder = new ProteinFinder(this.outDir, this.roleFile);
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the SOUR protein FASTA.
        log.info("Creating SOUR protein database.");
        this.finder.createProteinFile(this.genomes);
        // Create the reference-genome finder.
        log.info("Computing reference-genome FASTA files.");
        this.finder.createDnaFiles(this.refMap);
    }

}
