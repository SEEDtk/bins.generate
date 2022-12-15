/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.bins.BinGroup;
import org.theseed.bins.BinParms;

/**
 * This object represents a binning method.  It operates on a BinGroup with pre-selected starter bins, then
 * each method type uses its own rules to assign the leftover contigs to bins.
 *
 *
 * @author Bruce Parrello
 *
 */
public abstract class BinningMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinningMethod.class);

    /**
     * This enum contains the different binning method types.
     */
    public static enum Type {
        /** protein-based kmers for bin assignment */
        STANDARD {
            @Override
            public BinningMethod create(BinParms parms, File outDir) {
                // TODO code for create standard binning
                return null;
            }
        },
        /** do not bin, merely estimate species population */
        REPORT {
            @Override
            public BinningMethod create(BinParms parms, File outDir) {
                // TODO code for create report-only binning
                return null;
            }
        },
        /** DNA-based kmers for bin assignment */
        STRICT {
            @Override
            public BinningMethod create(BinParms parms, File outDir) {
                // TODO code for create DNA-based binning
                return null;
            }
        };

        /**
         * Create a binning engine for this type.
         *
         * @param parms		binning parameter object
         * @param outDir 	output directory for the bin data
         */
        public abstract BinningMethod create(BinParms parms, File outDir);

    }

    /**
     * Run the binning process on a bin group and save the group to the
     * specified output file.
     *
     * @param binGroup		binning group to process
     */
    public void run(BinGroup binGroup) {
        this.runMethod(binGroup);
    }

    /**
     * Use this method to assign contigs to bins.
     *
     * @param binGroup		binning group containing the contigs and starter bins
     */
    protected abstract void runMethod(BinGroup binGroup);

}
