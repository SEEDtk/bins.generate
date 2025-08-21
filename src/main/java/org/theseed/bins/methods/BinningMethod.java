/**
 *
 */
package org.theseed.bins.methods;

import java.util.List;

import org.theseed.bins.Bin;
import org.theseed.sequence.DiscriminatingKmerDb;
import org.theseed.sequence.FeatureDiscriminatingKmerDb;
import org.theseed.sequence.ProteinDiscriminatingKmerDb;

/**
 * This object represents a binning method.  It operates on a BinGroup with pre-selected starter bins, then
 * each method type uses its own rules to assign the leftover contigs to bins.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BinningMethod {

    // FIELDS
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
                DiscriminatingKmerDb db = new ProteinDiscriminatingKmerDb(processor.getParms().getKProt());
                return new KmerBinningMethod(processor, db);
            }
        },
        /** do not bin, merely estimate species population */
        REPORT {
            @Override
            public BinningMethod create(BinPhase.IParms processor) {
                return new NullBinningMethod(processor);
            }
        },
        /** DNA-based kmers for bin assignment */
        STRICT {
            @Override
            public BinningMethod create(BinPhase.IParms processor) {
                DiscriminatingKmerDb db = new FeatureDiscriminatingKmerDb(processor.getParms().getKDna());
                return new KmerBinningMethod(processor, db);
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
     *
     * @throws Exception
     */
    public void run(List<Bin> starterBins) throws Exception {
        this.runMethod(starterBins);
    }

    /**
     * Use this method to assign contigs to bins.
     *
     * @throws Exception
     */
    protected abstract void runMethod(List<Bin> starterBins) throws Exception;

}
