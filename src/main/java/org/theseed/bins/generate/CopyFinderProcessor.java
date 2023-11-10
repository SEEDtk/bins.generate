/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;


import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.utils.BaseInputProcessor;

/**
 * This command builds a finder by taking a subset of the roles in an existing finder.  The command uses a list
 * of roles IDs in the standard input to pull out the role subset from the source finder's role definition file,
 * and this guides what is copied from the source finder to the new finder.
 *
 * The oositional parameters are the name of the source finder and the name of the new, output finder.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more detailed log messages
 * -i	input file containing role IDs (if not STDIN)
 * -c	index (1-based) or name of input column containing role IDs (default "1")
 *
 * --clear	clear the output directory before processing
 *
 * The
 * @author Bruce Parrello
 *
 */
public class CopyFinderProcessor extends BaseInputProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CopyFinderProcessor.class);
    /** desired role subset */
    private Set<String> newRoles;
    /** input (source) protein finder */
    private ProteinFinder sourceFinder;
    /** map of role IDs to FASTA files in the input finder */
    private Map<String, File> fastaMap;

    // COMMAND-LINE OPTIONS

    /** input column containing role IDs */
    @Option(name = "--col", aliases = { "-c" }, metaVar = "role_id", usage = "index (1-based) or name of input column containing the IDs of the roles to use")
    private String colName;

    /** if specified, the output directory will be cleared before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** name of the source finder directory */
    @Argument(index = 0, metaVar = "inDir", usage = "name of the source (input) finder directory")
    private File inDir;

    /** name of the target finder directory */
    @Argument(index = 1, metaVar = "outDir", usage = "name of the new (output) finder directory")
    private File outDir;

    @Override
    protected void setReaderDefaults() {
        this.colName = "1";
        this.clearFlag = false;
    }

    @Override
    protected void validateReaderParms() throws IOException, ParseFailureException {
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Load the role definitions.
        this.sourceFinder = new ProteinFinder(this.inDir);
        this.fastaMap = this.sourceFinder.getFastas();
        // Set up the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating new output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("New finder will be created in {}.", this.outDir);
    }

    @Override
    protected void validateReaderInput(TabbedLineReader reader) throws IOException {
        // Get the column containing the role IDs.
        int roleColIdx = reader.findField(this.colName);
        // Read and validate the roles.
        this.newRoles = new TreeSet<String>();
        for (var line : reader) {
            String roleId = line.get(roleColIdx);
            if (! this.fastaMap.containsKey(roleId))
                throw new IOException("Role ID \"" + roleId + "\" from input is not in the finder at " + this.inDir + ".");
            this.newRoles.add(roleId);
        }
        log.info("{} of {} roles will be transferred to {}.", this.newRoles.size(), this.fastaMap.size(),
                this.outDir);
    }

    @Override
    protected void runReader(TabbedLineReader reader) throws Exception {
        // Copy the role map.
        File inRoleFile = this.sourceFinder.getRoleFile();
        File outRoleFile = new File(this.outDir, ProteinFinder.ROLE_FILE_NAME);
        log.info("Copying {} roles from {} to {}.", this.newRoles.size(), inRoleFile, outRoleFile);
        try (LineReader inRoles = new LineReader(inRoleFile);
                PrintWriter outRoles = new PrintWriter(outRoleFile)) {
            int inCount = 0;
            int outCount = 0;
            for (var line : inRoles) {
                inCount++;
                String inRole = StringUtils.substringBefore(line, "\t");
                if (this.newRoles.contains(inRole)) {
                    outCount++;
                    outRoles.println(line);
                }
            }
            log.info("{} role definition lines read, {} written.", inCount, outCount);
        }
        // Copy the protein FASTA.
        File inProtFile = new File(this.inDir, ProteinFinder.PROTEIN_FILE_NAME);
        File outProtFile = new File(this.outDir, ProteinFinder.PROTEIN_FILE_NAME);
        log.info("Copying seed proteins from {} to {}.", inProtFile, outProtFile);
        try (FastaInputStream inProts = new FastaInputStream(inProtFile);
                FastaOutputStream outProts = new FastaOutputStream(outProtFile)) {
            int inCount = 0;
            int outCount = 0;
            for (Sequence seq : inProts) {
                inCount++;
                if (this.newRoles.contains(seq.getComment())) {
                    outCount++;
                    outProts.write(seq);
                }
            }
            log.info("{} protein sequences read, {} written.", inCount, outCount);
        }
        int inCount = 0;
        int outCount = 0;
        // Finally, copy the role FASTA files.
        for (var roleFastaEntry : this.fastaMap.entrySet()) {
            inCount++;
            if (this.newRoles.contains(roleFastaEntry.getKey())) {
                final File inFastaFile = roleFastaEntry.getValue();
                log.info("Copying DNA fasta file {} to {}.", inFastaFile, this.outDir);
                FileUtils.copyFileToDirectory(inFastaFile, this.outDir, true);
                outCount++;
            }
        }
        log.info("{} of {} role FASTA files copied to {}.", outCount, inCount, this.outDir);
    }

 }
