/**
 *
 */
package org.theseed.bins.methods;

import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.bins.Bin;

/**
 * This object represents a binning method.  It operates on a BinGroup with pre-selected starter bins, then
 * each method type uses its own rules to assign the leftover contigs to bins.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BinningMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinningMethod.class);
    /** controlling command processor */
    protected BinPhase.IParms processor;

    /**
     * This enum contains the different binning method types.
     */
    public static enum Type {
        /** protein-based kmers for bin assignment */
        STANDARD {
            @Override
            public BinningMethod create(BinPhase.IParms processor) {
                // TODO code for create standard binning
                return null;
            }
        },
        /** do not bin, merely estimate species population */
        REPORT {
            @Override
            public BinningMethod create(BinPhase.IParms processor) {
                // TODO code for create report-only binning
                return null;
            }
        },
        /** DNA-based kmers for bin assignment */
        STRICT {
            @Override
            public BinningMethod create(BinPhase.IParms processor) {
                // TODO code for create DNA-based binning
                return null;
            }
        };

        /**
         * Create a binning engine for this type.
         *
         * @processor	controlling command processor
         */
        public abstract BinningMethod create(BinPhase.IParms processor);

    }

    /**
     * Construct a binning method.
     *
     * @param controller		controlling command processor
     */
    public BinningMethod(BinPhase.IParms controller) {
        this.processor = controller;
    }

    /**
     * Run the binning process.
     *
     * @param starterBins	list of significant bins set up as starters
     */
    public void run(List<Bin> starterBins) {
        this.runMethod(starterBins);
    }

    /**
     * Use this method to assign contigs to bins.
     */
    protected abstract void runMethod(List<Bin> starterBins);

}
