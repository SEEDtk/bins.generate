package org.theseed.bins.generate;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;

/**
 * Commands for utilities relating to binning.
 *
 * build		build the protein finder FASTA files
 * copy			create a new protein finder that is a subset of a bigger one
 * bin			process a FASTA file to create bins
 * clean		remove ambiguous sequences from a finder's FASTA files
 * sourFile		create a subset of a role definition file with a specified set of roles
 * checkv_db	update the checkv database to include taxon IDs
 *
 */
public class App {

    /** static array containing command names and comments */
    protected static final String[] COMMANDS = new String[] {
             "build", "build the protein finder FASTA files",
             "copy", "create a new protein finder that is a subset of a bigger one",
             "bin", "process a FASTA file to create bins",
             "clean", "remove ambiguous sequences from a finder's FASTA files",
             "sourFile", "create a subset of a role definition file with a specified set of roles",
             "checkv_db", "update the checkv database to include taxon IDs",
    };

    public static void main( String[] args ) {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Determine the command to process.
        switch (command) {
        case "build" :
            processor = new BuildFinderProcessor();
            break;
        case "copy" :
            processor = new CopyFinderProcessor();
            break;
        case "bin" :
            processor = new BinProcessor();
            break;
        case "clean" :
            processor = new CleanFinderProcessor();
            break;
        case "sourFile" :
            processor = new SourFileProcessor();
            break;
        case "checkv_db" :
            processor = new CheckVDbProcessor();
            break;
        case "-h" :
        case "--help" :
            processor = null;
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        if (processor == null)
            BaseProcessor.showCommands(COMMANDS);
        else {
            processor.parseCommand(newArgs);
            processor.run();
        }
    }
}
