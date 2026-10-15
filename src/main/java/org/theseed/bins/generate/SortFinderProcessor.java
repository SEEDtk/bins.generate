package org.theseed.bins.generate;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.kohsuke.args4j.Argument;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This command sorts the FASTA files in a protein finder. Each FASTA is loaded into memory (a typical finder FASTA is 2G) 
 * and then written to a temp file sorted by label. The temp file then replaces the original FASTA file. The resulting
 * finders can then be used to create a FinderKmerBatch save file for the entire good-genome set.
 * 
 * The positional parameter is the name of the protein finder directory. The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 */
public class SortFinderProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(SortFinderProcessor.class);
    /** protein finder manager */
    private ProteinFinder finder;

    // COMMAND-LINE OPTIONS

    /** protein finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory", required = true)
    private File finderDir;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        if (! this.finderDir.isDirectory())
            throw new IOException("Finder directory " + this.finderDir + " is not found or invalid.");
        this.finder = new ProteinFinder(this.finderDir);
    }
    
    @Override
    protected void runCommand() throws Exception {
        // Get the map of role IDs to FASTA files.
        Map<String, File> roleMap = this.finder.getFastas();
        final int nFiles = roleMap.size();
        log.info("{} FASTA files found.", nFiles);
        // This will count the number of files processed.
        int processed = 0;
        // Loop through the files, sorting them.
        for (File fastaFile : roleMap.values()) {
            processed++;
            log.info("Processing file {} of {}: {}", processed, nFiles, fastaFile);
            // Load the FASTA file into memory. The tree map will automatically sort the sequences by their labels.
            TreeMap<String, Sequence> seqMap = new TreeMap<>();
            try (FastaInputStream fastaIn = new FastaInputStream(fastaFile)) {
                fastaIn.forEach(seq -> seqMap.put(seq.getLabel(), seq));
                log.info("Loaded {} sequences from {}", seqMap.size(), fastaFile);
            }
            // Write the sorted sequences to a temporary file.
            File tempFile = new File(fastaFile.getParentFile(), fastaFile.getName() + ".tmp");
            try (FastaOutputStream fastaOut = new FastaOutputStream(tempFile)) {
                for (Sequence seq : seqMap.values()) {
                    fastaOut.write(seq);
                }
            }
            log.info("Renaming temporary file {} to {}", tempFile, fastaFile);
            if (! fastaFile.delete())
                throw new IOException("Could not delete original FASTA file " + fastaFile);
            if (! tempFile.renameTo(fastaFile))
                throw new IOException("Could not rename temp file " + tempFile + " to " + fastaFile);
        }
        log.info("{} FASTA files sorted successfully.", nFiles);
    }
    
}
