/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;

/**
 * This command takes as input a master role definition file and produces a role definition file for a selected subset of roles.
 *
 * The positional parameter is the name of a tab-delimited file (with headers) containing the IDs of the roles to keep.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input role definition file (if not STDIN)
 * -o 	output role definition file (if not STDOUT)
 *
 * --col	index (1-based) or name of the role-set-file column containing role IDs (default "1")
 *
 * @author Bruce Parrello
 *
 */
public class SourFileProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SourFileProcessor.class);
    /** set of role IDs to keep */
    private Set<String> keepIdSet;

    // COMMAND-LINE OPTIONS

    /** name of the input role definition file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "roles.in.subsystems", usage = "name of the input role definition file (if not STDIN)")
    private File roleFileName;

    /** index (1-based) or name of the role-set-file column containing roles IDs */
    @Option(name = "--col", metaVar = "role", usage = "index (1-based) or name of the role-set-file column containing roles IDs")
    private String colName;

    /** role-set file name */
    @Argument(index = 0, metaVar = "roleSetFile", usage = "name of the file containing the IDs of the roles to keep", required = true)
    private File roleSetFile;

    @Override
    protected void setReporterDefaults() {
        this.colName = "1";
        this.roleFileName = null;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that the input file is readable if one was specified.
        if (this.roleFileName != null && ! this.roleFileName.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFileName + " is not found or unreadable.");
        // We need to read in the role IDs to keep.  Verify that the role-set file exists.
        if (! this.roleSetFile.canRead())
            throw new FileNotFoundException("Role set file " + this.roleSetFile + " is not found or unreadable.");
        log.info("Reading role IDs from column \"{}\" of {}.", this.colName, this.roleSetFile);
        this.keepIdSet = TabbedLineReader.readSet(this.roleSetFile, this.colName);
        log.info("{} roles will be kept.", this.keepIdSet.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Open the input file.
        try (LineReader inStream = new LineReader(this.roleFileName)) {
            log.info("Input role definitions taken from {}.", inStream.getName());
            int inCount = 0;
            int outCount = 0;
            // Loop through the role definitions.
            for (String line : inStream) {
                inCount++;
                String roleId = StringUtils.substringBefore(line, "\t");
                if (this.keepIdSet.contains(roleId)) {
                    outCount++;
                    writer.println(line);
                }
            }
            log.info("{} lines read.  {} lines written.", inCount, outCount);
        }
    }


}
