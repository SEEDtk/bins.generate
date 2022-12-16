/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.theseed.bins.Bin;
import org.theseed.bins.BinGroup;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.sequence.DiscriminatingKmerDb;
import org.theseed.sequence.DnaDiscriminatingKmerDb;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This is the basic binning method, which uses discriminating kmers.  The type of kmer database is
 * determined by the client.  Each reference genome is added to the database and then discriminating
 * kmers are counted in each unplaced contig to determine which bin it should be put in.
 *
 * @author Bruce Parrello
 *
 */
public class KmerBinningMethod extends BinningMethod {

    // FIELDS
    /** discriminating-kmer database */
    private DiscriminatingKmerDb kmerDb;
    /** map of bin names to bins */
    private Map<String, Bin> binMap;

    /**
     * Construct a kmer-based binning object.  The database is built and destroyed by this object.
     *
     * @param processor		controlling command processor
     * @param db			kmer database to use
     */
    public KmerBinningMethod(BinPhase.IParms processor, DiscriminatingKmerDb db) {
        super(processor);
        this.kmerDb = db;
    }

    @Override
    protected void runMethod(List<Bin> starterBins) throws IOException {
        // Get a map of bin IDs to starter bins.  We will use the bin IDs as group IDs in the kmer map.
        this.binMap = starterBins.stream().collect(Collectors.toMap(x -> x.getID(), x -> x));
        // Make an initial pass using reference genomes.
        this.processRefGenomes();
        // Check the repeat-region parameter.
        int dangLen = this.processor.getParms().getDangLen();
        if (dangLen > 0) {
            // Make a second pass using long contig kmers from the bins.
            this.processRepeatRegions(dangLen);
        }
    }

    /**
     * This method builds the main kmer database from the reference genomes in the starter bins, and then
     * uses kmer counts to assign contigs to bins.
     *
     * @throws IOException
     */
    private void processRefGenomes() throws IOException {
        // Loop through the starter bins.
        for (Bin bin : this.binMap.values()) {
            // The group name is the name of the first contig in the bin.  This allows us to find
            // the bin using the contig map.
            String name = bin.getContigs().get(0);
            log.info("Processing reference genomes for starter bin {} using ID {}.", bin.getName(), name);
            // Get the reference genomes and extract the kmers.
            Collection<String> refs = bin.getAllRefGenomes();
            for (String refGenomeId : refs) {
                Genome refGenome = processor.getGenome(refGenomeId);
                log.info("Scanning for kmers in {}.", refGenome);
                this.kmerDb.addGenome(refGenome, name);
            }
        }
        // Finalize the kmers.
        this.kmerDb.finalize();
        // Get the binning group.
        BinGroup binGroup = this.processor.getBinGroup();
        // Get the bin strength limit.  This is the minimum difference required from the best hit count and
        // the second-best hit count.
        int minDifference = this.processor.getParms().getBinStrength();
        // There are a few counts we want to display for tracing.
        int contigCount = 0;
        int placeCount = 0;
        // Now we read the contigs and assign them to bins.
        File inFile = this.processor.getInFile();
        log.info("Scanning contigs in {}.", inFile);
        try (FastaInputStream inStream = new FastaInputStream(inFile)) {
            for (Sequence contig : inStream) {
                contigCount++;
                // Only proceed if this contig has a bin and is unplaced.
                Bin contigBin = binGroup.getContigBin(contig.getLabel());
                if (contigBin == null) {
                    // Here the bin was rejected because of low quality.
                    binGroup.count("kmer-contig-Skip");
                } else if (! contigBin.isSignificant()) {
                    // Here the bin is of sufficient quality and is not in a bin yet.  Scan it for kmers.
                    CountMap<String> counts = this.kmerDb.countHits(contig.getSequence());
                    // Get the counts in order.
                    var counters = counts.sortedCounts();
                    if (counters.size() <= 0) {
                        // No kmer hits were found.  We have no idea what this is.
                        binGroup.count("kmer-contig-NoHits");
                    } else {
                        String targetBinName = counters.get(0).getKey();
                        int targetCount = counters.get(0).getCount();
                        if (counters.size() < 2 ) {
                            // Only one group works for this contig.
                            if (targetCount < minDifference) {
                                // Not enough hits, so treat it as noise.
                                binGroup.count("kmer-contig-UnambiguousWeak");
                                targetBinName = null;
                            } else
                                binGroup.count("kmer-contig-UnambiguousStrong");
                        } else {
                            // Here we have multiple counters.  Make sure the difference is high enough.
                            if (targetCount - counters.get(1).getCount() < minDifference) {
                                // The counts are too close.  Reject this contig.
                                binGroup.count("kmer-contig-AmbiguousWeak");
                                targetBinName = null;
                            } else
                                binGroup.count("kmer-contig-AmbiguousStrong");
                        }
                        if (targetBinName != null) {
                            // Merge the new contig bin into the target bin.
                            Bin targetBin = binGroup.getContigBin(targetBinName);
                            binGroup.merge(targetBin, contigBin);
                            binGroup.count("kmer-contig-Placed");
                            placeCount++;
                        }
                    }
                }
                if (log.isInfoEnabled() && contigCount % 1000 == 0) {
                    log.info("{} contigs read, {} placed.", contigCount, placeCount);
                }
            }
        }
        log.info("{} contigs read and {} placed by discriminating-kmer analysis.", contigCount, placeCount);
        // Erase the discriminating-kmer database to save memory.
        this.kmerDb.clear();
    }

    /**
     * This method extracts long DNA sequences from the contigs already in bins, and looks for them in the
     * unplaced contigs.  This requires a staggering two passes over the input FASTA file.
     *
     * @param dangLen	kmer size to use
     *
     * @throws IOException
     */
    private void processRepeatRegions(int dangLen) throws IOException {
        // Get some counters for tracing.
        int readCount = 0;
        int scanCount = 0;
        int placeCount = 0;
        // Get the binning group.
        BinGroup binGroup = this.processor.getBinGroup();
        // Build a contig-based DNA kmer database and populate it with kmers from the bins.
        DnaDiscriminatingKmerDb dangKmers = new DnaDiscriminatingKmerDb(dangLen);
        File inFile = this.processor.getInFile();
        log.info("Scanning {} for repeat-region kmers.", inFile);
        try (FastaInputStream inStream = new FastaInputStream(inFile)) {
            for (Sequence seq : inStream) {
                readCount++;
                Bin contigBin = binGroup.getContigBin(seq.getLabel());
                if (contigBin != null && contigBin.isSignificant()) {
                    // This contig is in a bin, so scan it.
                    scanCount++;
                    dangKmers.addContig(contigBin.getName(), seq.getSequence());
                }
                if (log.isInfoEnabled() && readCount % 1000 == 0)
                    log.info("{} contigs read and {} scanned for repeat-region kmers.", readCount, scanCount);
            }
        }
        // Finalize the database.
        dangKmers.finalize();
        // Now we make the second pass, looking for hits in the unplaced contigs.
        readCount = 0;
        scanCount = 0;
        log.info("Scanning {} for repeat-region placement.", inFile);
        try (FastaInputStream inStream = new FastaInputStream(inFile)) {
            for (Sequence seq : inStream) {
                readCount++;
                Bin contigBin = binGroup.getContigBin(seq.getLabel());
                if (contigBin != null && ! contigBin.isSignificant()) {
                    // This contig is not in a bin, so check it.
                    scanCount++;
                    CountMap<String> counts = dangKmers.countHits(seq.getSequence());
                    var bestCount = counts.getBestEntry();
                    if (bestCount != null) {
                        // Here we found a likely repeat region.  Figure out where to place it.
                        Bin targetBin = this.binMap.get(bestCount.getKey());
                        binGroup.merge(targetBin, contigBin);
                        placeCount++;
                        binGroup.count("repeat-contig-Placed");
                    } else
                        binGroup.count("repeat-contig-NoHits");
                }
                if (log.isInfoEnabled() && readCount % 1000 == 0)
                    log.info("{} contigs read, {} checked, and {} placed during repeat-region placement pass.",
                            readCount, scanCount, placeCount);
            }
        }
    }

}
