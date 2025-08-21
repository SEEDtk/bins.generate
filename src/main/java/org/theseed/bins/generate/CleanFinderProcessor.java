/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.sequence.seeds.FastaCleaner;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This is a repair method for finders.  We copy each FASTA file, removing ambiguous sequences.
 *
 * The positional parameter is the name of the finder directory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * @author Bruce Parrello
 *
 */
public class CleanFinderProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(CleanFinderProcessor.class);
    /** protein finder */
    private ProteinFinder finder;

    // COMMAND-LINE OPTIONS

    /** name of the finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "name of the finder directory")
    private File inDir;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        // Insure the finder directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input finder directory " + this.inDir + " is not found or invalid.");
        log.info("Loading protein finder from {}.", this.inDir);
        this.finder = new ProteinFinder(this.inDir);
    }

    @Override
    protected void runCommand() throws Exception {
        // First, clean the protein file.
        File protFile = this.finder.getProteinFile();
        FastaCleaner cleaner = new FastaCleaner.Protein();
        log.info("Processing protein file {}.", protFile);
        cleaner.clean(protFile);
        // Create a DNA cleaner.
        cleaner = new FastaCleaner.Dna();
        // Now clean the DNA files.
        Map<String, File> roleFiles = this.finder.getFastas();
        for (var roleEntry : roleFiles.entrySet()) {
            String roleId = roleEntry.getKey();
            File roleFile = roleEntry.getValue();
            // If the file exists, clean it. If it is not found, build it.
            if (roleFile.exists()) {
                log.info("Processing DNA file {} for {}.", roleFile, roleId);
                cleaner.clean(roleFile);
            }
        }
    }

}
