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
    /** minimum length for a SOUR protein search contig */
    private int lenFilter;
    /** minimum length for a binning contig */
    private int binLenFilter;
    /** minimum coverage for a SOUR protein search contig */
    private double covgFilter;
    /** minimum coverage for a binning contig */
    private double binCovgFilter;
    /** maximum number of ambiguity characters in a row for a contig */
    private int xLimit;
    /** maximum e-value for BLAST hits during the SOUR protein search */
    private double maxEValue;
    /** maximum e-value for BLAST hits during the reference-genome search */
    private double refMaxEValue;
    /** minimum fraction of a protein that must match in a BLAST hit */
    private double minLen;
    /** maximum allowable gap for merging during the SOUR protein search */
    private int maxGap;
    /** protein kmer length */
    private int kProt;
    // DNA kmer length */
    private int kDna;
    /** repeat-region kmer length */
    private int dangLen;
    /** hit-differential for contig assignment */
    private int binStrength;

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
        this.setMaxGap(600);
        this.setKProt(8);
        this.setKDna(15);
        this.setDangLen(50);
        this.setBinStrength(10);
    }

    /**
     * @return the maximum gap for hit-merging in the SOUR protein search
     */
    public int getMaxGap() {
        return this.maxGap;
    }

    /**
     * Specify a new maximum gap for hit-merging in the SOUR protein search.
     *
     * @param maxGap	the maxGap to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setMaxGap(int maxGap) {
        this.maxGap = maxGap;
        return this;
    }

    /**
     * @return the minimum length for a SOUR protein search contig
     */
    public int getLenFilter() {
        return this.lenFilter;
    }

    /**
     * Specify the minimum length for a SOUR protein search contig.
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
     * @return the minimum coverage for a SOUR protein search contig
     */
    public double getCovgFilter() {
        return this.covgFilter;
    }

    /**
     * Specify the minimum coverage for a SOUR protein search contig
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
     * @return the maximum e-value for BLAST hits during the SOUR protein search
     */
    public double getMaxEValue() {
        return this.maxEValue;
    }

    /**
     * Specify the maximum e-value for BLAST hits during the SOUR protein search.
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

    /**
     * @return the protein kmer size
     */
    public int getKProt() {
        return this.kProt;
    }

    /**
     * Specify a new protein kmer size.
     *
     * @param kProt 	the kProt to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setKProt(int kProt) {
        this.kProt = kProt;
        return this;
    }

    /**
     * @return the DNA kmer size
     */
    public int getKDna() {
        return this.kDna;
    }

    /**
     * Specify a new DNA kmer size.
     *
     * @param kDna the kDna to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setKDna(int kDna) {
        this.kDna = kDna;
        return this;
    }

    /**
     * @return the repeat-region kmer length
     */
    public int getDangLen() {
        return this.dangLen;
    }

    /**
     * Specify a new repeat-region kmer length.
     *
     * @param danglen the dangLen to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setDangLen(int danglen) {
        this.dangLen = danglen;
        return this;
    }

    /**
     * @return the hit-differential minimum
     */
    public int getBinStrength() {
        return this.binStrength;
    }

    /**
     * Specify a new hit-differential minimum.
     *
     * @param binStrength	the bin strength to set
     *
     * @return this object, for fluent invocation
     */
    public BinParms setBinStrength(int binStrength) {
        this.binStrength = binStrength;
        return this;
    }


    @Override
    public String toString() {
        return String.format(
                "--lenFilter=%s --binLenFilter=%s --covgFilter=%s --binCovgFilter=%s, --xLimit=%s  --maxEValue=%s" +
                " --refMaxEValue=%s --minLen=%s --maxGap=%s --kProt=%s --kDna=%s --dangLen=%s --binStrength=%s",
                this.lenFilter, this.binLenFilter, this.covgFilter, this.binCovgFilter, this.xLimit, this.maxEValue,
                this.refMaxEValue, this.minLen, this.maxGap, this.kProt, this.kDna, this.dangLen, this.binStrength);
    }

}
