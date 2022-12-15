/**
 *
 */
package org.theseed.bins.methods;

import java.util.List;

import org.theseed.bins.Bin;
import org.theseed.bins.methods.BinPhase.IParms;

/**
 * This binning method does not assign any contigs to bins, so the output report only contains the starter bins.
 *
 * @author Bruce Parrello
 *
 */
public class NullBinningMethod extends BinningMethod {

    /**
     * Construct a null binner.
     *
     * @param controller	controlling command processor
     */
    public NullBinningMethod(IParms controller) {
        super(controller);
    }

    @Override
    protected void runMethod(List<Bin> starterBins) {
    }

}
