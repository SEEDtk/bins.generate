/**
 *
 */
package org.theseed.bins.methods;

/**
 * This phase actually assigns contigs to the significant bins.  Most of the work is done by the
 * binning engine.
 *
 * @author Bruce Parrello
 *
 */
public class RunPhase extends BinPhase {

    public RunPhase(IParms commandProcessor) {
        super(commandProcessor);
    }

    @Override
    protected String getSaveFileName() {
        return "bins.kmers.json";
    }

    @Override
    protected String getPhaseName() {
        return "CONTIG-ASSIGNMENT";
    }

    @Override
    protected void runPhase() {
        // Get the significant bins.
        var starterBins = this.getBinGroup().getSignificantBins();
        // Run the binning method if there are any.
        if (! starterBins.isEmpty())
            this.getEngine().run(starterBins);
    }

}
