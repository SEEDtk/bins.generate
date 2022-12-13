/**
 *
 */
package org.theseed.bins;

/**
 * This class manages the basic parameters for binning.
 *
 * @author Bruce Parrello
 */
public class BinParms {

    // FIELDS
    /** minimum length for a seed-protein search contig */
    private int lenFilter;
    /** minimum length for a binning contig */
    private int binLenFilter;
    /** minimum coverage for a seed-protein search contig */
    private double covgFilter;
    /** minimum coverage for a binning contig */
    private double binCovgFilter;
    /** maximum number of ambiguity characters in a row for a contig */
    private int xLimit;
    /** maximum e-value for BLAST hits during the seed-protein search */
    private double maxEValue;
    /** maximum e-value for BLAST hits during the reference-genome search */
    private double refMaxEValue;
    /** minimum fraction of a protein that must match in a BLAST hit */
    private double minLen;


    /**
     * Generate a binning parameter set with default parameters.
     */
    public BinParms() {
        this.setBinCovgFilter(5.0);
        this.setBinLenFilter(300);
        this.setCovgFilter(5.0);
        this.setLenFilter(500);
        this.setXLimit(10);
        this.setMaxEValue(1e-20);
        this.setRefMaxEValue(1e-10);
        this.setMinLen(0.5);
    }

    /**
     * @return the minimum length for a seed-protein search contig
     */
    public int getLenFilter() {
        return this.lenFilter;
    }

    /**
     * Specify the minimum length for a seed-protein search contig.
     *
     * @param lenFilter the lenFilter to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setLenFilter(int lenFilter) {
        this.lenFilter = lenFilter;
        return this;
    }

    /**
     * @return the minimum length for a binning contig
     */
    public int getBinLenFilter() {
        return this.binLenFilter;
    }

    /**
     * Specify the minimum length for a binning contig
     *
     * @param binLenFilter the binLenFilter to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setBinLenFilter(int binLenFilter) {
        this.binLenFilter = binLenFilter;
        return this;
    }

    /**
     * @return the minimum coverage for a seed-protein search contig
     */
    public double getCovgFilter() {
        return this.covgFilter;
    }

    /**
     * Specify the minimum coverage for a seed-protein search contig
     *
     * @param covgFilter the covgFilter to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setCovgFilter(double covgFilter) {
        this.covgFilter = covgFilter;
        return this;
    }

    /**
     * @return the minimum coverage for a binning contig
     */
    public double getBinCovgFilter() {
        return this.binCovgFilter;
    }

    /**
     * Specify the minimum coverage for a binning contig.
     *
     * @param binCovgFilter the binCovgFilter to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setBinCovgFilter(double binCovgFilter) {
        this.binCovgFilter = binCovgFilter;
        return this;
    }

    /**
     * @return the maximum number of ambiguity characters in a row for a contig
     */
    public int getXLimit() {
        return this.xLimit;
    }

    /**
     * Specify the maximum number of ambiguity characters in a row for a contig.
     *
     * @param xLimit the xLimit to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setXLimit(int xLimit) {
        this.xLimit = xLimit;
        return this;
    }

    /**
     * @return the maximum e-value for BLAST hits during the seed-protein search
     */
    public double getMaxEValue() {
        return this.maxEValue;
    }

    /**
     * Specify the maximum e-value for BLAST hits during the seed-protein search.
     *
     * @param maxEValue the maxEValue to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setMaxEValue(double maxEValue) {
        this.maxEValue = maxEValue;
        return this;
    }

    /**
     * @return the maximum e-value for BLAST hits during the reference-genome search
     */
    public double getRefMaxEValue() {
        return this.refMaxEValue;
    }

    /**
     * Specify the maximum e-value for BLAST hits during the reference-genome search
     *
     * @param refMaxEValue the refMaxEValue to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setRefMaxEValue(double refMaxEValue) {
        this.refMaxEValue = refMaxEValue;
        return this;
    }

    /**
     * @return the minimum fraction of a protein that must match in a BLAST hit
     */
    public double getMinLen() {
        return this.minLen;
    }

    /**
     * Specify the minimum fraction of a protein that must match in a BLAST hit.
     *
     * @param minLen the minLen to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setMinLen(double minLen) {
        this.minLen = minLen;
        return this;
    }


}
