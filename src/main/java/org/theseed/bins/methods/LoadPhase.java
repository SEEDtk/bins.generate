/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;
import java.io.IOException;

import org.theseed.bins.BinGroup;

/**
 * This is the first phase of binning.  It loads the sample's FASTA file into the binning
 * group and creates the reduced FASTA file used for the SOUR search.
 *
 * @author Bruce Parrello
 *
 */
public class LoadPhase extends BinPhase {

    // FIELDS
    /** name of the reduced FASTA file used for the SOUR search */
    public static final String REDUCED_FASTA_NAME = "reduced.fasta";

    public LoadPhase(IParms commandProcessor) {
        super(commandProcessor);
    }

    @Override
    protected String getSaveFileName() {
        return "bin.contigs.json";
    }

    @Override
    protected String getPhaseName() {
        return "CONTIG-LOAD";
    }

    @Override
    protected void runPhase() throws IOException {
        // Get the name of the reduced-fasta file.  We put the SOUR protein search contigs in here.
        File reducedFile = this.getOutFile(REDUCED_FASTA_NAME);
        // Load the contigs into the binning group.
        BinGroup binGroup = this.getBinGroup();
        binGroup.loadFromFasta(this.getInFile(), this.getParms(), reducedFile);

    }

}
