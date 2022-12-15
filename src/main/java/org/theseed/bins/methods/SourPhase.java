/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.collections4.iterators.ReverseListIterator;
import org.theseed.bins.Bin;
import org.theseed.bins.BinGroup;
import org.theseed.genome.Genome;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This binning phase uses a protein finder to select starter bins for each identifiable species in the sample.  The
 * reduced fasta file produced by the load phase is scanned for likely SOUR proteins, and then these are associated
 * with reference genomes.  Contigs belonging to the same species are grouped together and the best reference genome
 * selected.
 *
 * @author Bruce Parrello
 *
 */
public class SourPhase extends BinPhase {

    // FIELDS
    /** name of the save file for the sour hit map */
    private static final String SOUR_MAP_NAME = "sours.found.tbl";
    /** name of the save file for the reference genome map */
    private static final String REF_GENOME_MAP_NAME = "ref.genomes.tbl";
    /** sorter for reference-genome hits by score */
    private static final Comparator<ProteinFinder.DnaHit> SCORE_SORTER = new ProteinFinder.DnaHit.ScoreSorter();

    /**
     * @param commandProcessor
     */
    public SourPhase(IParms commandProcessor) {
        super(commandProcessor);
    }

    @Override
    protected String getSaveFileName() {
        return "bins.starter.json";
    }

    @Override
    protected String getPhaseName() {
        return "SOUR-SEARCH";
    }

    @Override
    protected void runPhase() throws IOException, InterruptedException {
        // Get the bin group.
        BinGroup binGroup = this.getBinGroup();
        // Find the SOUR proteins in the reduced FASTA file output by the load phase.
        File reducedFastaFile = this.getOutFile(LoadPhase.REDUCED_FASTA_NAME);
        ProteinFinder finder = this.getFinder();
        log.info("Searching for SOUR proteins in {}.", reducedFastaFile);
        var soursFound = finder.findSourProteins(reducedFastaFile);
        if (soursFound.size() <= 0) {
            // Here no SOUR proteins were found.
            log.warn("No SOUR proteins could be found in this sample.");
        } else {
            binGroup.count("sour-proteins-found", soursFound.size());
            // Checkpoint the results for debugging purposes.
            ProteinFinder.saveSeedProteins(soursFound, this.getOutFile(SOUR_MAP_NAME));
            // Now assign the reference genomes.
            log.info("Searching for reference genomes for {} SOUR protein regions.", soursFound.size());
            var refsFound = finder.findRefGenomes(soursFound, reducedFastaFile);
            if (refsFound.size() <= 0) {
                // Here we can't bin, because everything is too far away from our reference genome set.
                log.warn("No reference genomes could be found for this sample.");
            } else {
                binGroup.count("sour-refGenomes-assigned", refsFound.size());
                // Checkpoint these results as well.
                ProteinFinder.saveRefGenomes(refsFound, this.getOutFile(REF_GENOME_MAP_NAME));
                // Now we have a reference genome assigned to each contig, and it is associated with a score.
                // Both of these are stored in the DnaHit object.  We need to organize the DnaHits by species.
                var speciesMap = new HashMap<Integer, List<ProteinFinder.DnaHit>>(refsFound.size() * 4 / 3 + 1);
                for (ProteinFinder.DnaHit refHit : refsFound.values()) {
                    int taxId = refHit.getSpeciesId();
                    List<ProteinFinder.DnaHit> hitList = speciesMap.computeIfAbsent(taxId, x -> new ArrayList<ProteinFinder.DnaHit>(5));
                    hitList.add(refHit);
                    binGroup.count("sour-hit-" + refHit.getRoleId());
                }
                // For each species we need to pick the best reference genome and assign it to all the contigs, which are then merged
                // into a single per-species bin.  We need the name suffix so we can name the bins.
                String nameSuffix = this.getNameSuffix();
                // Loop through the species lists found.
                for (var speciesHitList : speciesMap.values()) {
                    Collections.sort(speciesHitList, SCORE_SORTER);
                    // We want to iterate through this list from the best-scoring hit to the worst, so we use a ReversListIterator.
                    var iter = new ReverseListIterator<ProteinFinder.DnaHit>(speciesHitList);
                    // Start with the best hit.  This determines the reference genome and becomes the master bin.
                    var bestHit = iter.next();
                    String refGenomeId = bestHit.getRefId();
                    Genome refGenome = this.getGenome(refGenomeId);
                    Bin masterBin = binGroup.getContigBin(bestHit.getContigId());
                    masterBin.setTaxInfo(bestHit, nameSuffix, refGenome);
                    binGroup.count("sour-species-found");
                    // Loop through the remaining hits, merging the contigs into the master.
                    while (iter.hasNext()) {
                        var hit = iter.next();
                        Bin bin = binGroup.getContigBin(hit.getContigId());
                        binGroup.merge(masterBin, bin);
                        masterBin.addRefGenome(hit.getRefId());
                        binGroup.count("sour-contig-merged");
                    }
                    log.info("{} contigs in starter bin {} using reference genome {}.", speciesHitList.size(), bestHit.getContigId(),
                            refGenome);
                }
            }
        }
    }

}
