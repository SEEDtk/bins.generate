/**
 *
 */
package org.theseed.bins.generate;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.bins.BinGroup;
import org.theseed.bins.BinParms;
import org.theseed.p3api.P3Connection;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 *
 * This is the new binning command.  Binning involves four main phases.
 *
 * 	1.	Filtering the contigs and computing coverage to build a BinGroup.
 *  2.	Locating seed proteins and building the starter bins.
 *  3.	Assigning the remaining contigs to bins.
 *  4.	Writing the individual bins to FASTA files.
 *
 *  We allow multiple possible methods for the third phase.  A full binning pipeline includes assembly of reads into
 *  contigs and annotating and evaluating the bins, but this command only does the middle part, which is the binning itself.
 *
 * The positional parameters are the name of the input FASTA file, and the name of the output directory for the binning.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * The following command-line options relate to the determination of the initial bins.
 *
 * --finder			name of the directory containing the protein finder files (default "Finder" in the current directory
 * --lenFilter		minimum length of a contig to be considered for the seed-protein search (default 300)
 * --binCovgFilter	minimum coverage for a contig to be considered for the seed-protein search (default 5.0)
 * --maxE			maximum e-value for BLAST hits when finding a seed protein (default 1e-20)
 * --refMaxE		maximum e-value for BLAST hits when finding a reference genome using a seed protein (default 1e-10)
 * --minLen			minimum fraction of a protein that must match in a protein BLAST hit (default 0.5)
 *
 * The following command-line options affect the entire binning process.
 *
 * --binLenFilter	minimum length of a contig to be considered for binning (default 300)
 * --binCovgFilter	minimum coverage for a contig to be considered for binning (default 5.0)
 * --xLimit			maximum run length of ambiguity characters allowed in a contig to be considered for binning
 * 					(default 10)
 * --kProt			protein kmer size (default 8)
 * --kDna			DNA kmer size (default 15)
 * --danglen		mobile element kmer size (default 50)
 * --recipe			method to be used for binning (default STANDARD)
 *
 * The following command-line options relate to the PATRIC database.
 *
 * --dataAPIUrl		alternate URL for accessing the PATRIC data service used to download reference genomes
 * 					(default from environment variable "P3API_URL")
 * --suffix			suffix to append to species names when computing bin names (default "clonal population")
 *
 * @author Bruce Parrello
 *
 */
public class BinProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinProcessor.class);
    /** binning parameter object */
    private BinParms parms;
    /** connection to PATRIC */
    private P3Connection p3;
    /** master bin group */
    private BinGroup binGroup;
    /** protein finder */
    private ProteinFinder finder;

    @Override
    protected void setDefaults() {
        // TODO code for setDefaults

    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // TODO code for validateParms
        return false;
    }

    @Override
    protected void runCommand() throws Exception {
        // TODO code for runCommand

    }
    // FIELDS
    // TODO data members for BinProcessor

    // TODO constructors and methods for BinProcessor
}
