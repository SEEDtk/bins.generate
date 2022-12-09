package org.theseed.bins.generate;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * Commands for utilities relating to binning.
 *
 * build		build the protein finder FASTA files
 *
 */
public class App {

    public static void main( String[] args ) {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Determine the command to process.
        switch (command) {
        case "build" :
            processor = new BuildProcessor();
            break;
        default:
            throw new RuntimeException("Invalid command " + command);
        }
        // Process it.
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
