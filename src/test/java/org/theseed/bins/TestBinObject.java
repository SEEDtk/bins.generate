/**
 *
 */
package org.theseed.bins;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;

import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * @author Bruce Parrello
 *
 */
class TestBinObject {

    @Test
    void testBinManipulation() throws IOException, JsonException {
        // Read in some useful DNA hits.
        var hits = ProteinFinder.loadRefGenomes(new File("data", "hitList.tbl"));
        // Create a bin group.
        BinGroup binGroup = new BinGroup();
        // Create a simple bin.
        final String node1Id = "NODE_58_length_26115_cov_109.467_ID_1122748";
        Bin bin1 = new Bin(node1Id, 1122748, 109.467);
        assertThat(bin1.isSignificant(), equalTo(false));
        assertThat(bin1.getContigs(), contains(node1Id));
        assertThat(bin1.getCoverage(), closeTo(109.467, 0.001));
        assertThat(bin1.getLen(), equalTo(1122748));
        // Make it significant.
        var hit = hits.get(node1Id);
        Genome refGenome = new Genome(new File("data", "548.1022.gto"));
        bin1.setTaxInfo(hit, "test bin", refGenome);
        assertThat(bin1.isSignificant(), equalTo(true));
        assertThat(bin1.getContigs(), contains(node1Id));
        assertThat(bin1.getCoverage(), closeTo(109.467, 0.001));
        assertThat(bin1.getLen(), equalTo(1122748));
        assertThat(bin1.getTaxonID(), equalTo(548));
        assertThat(bin1.getName(), equalTo("Klebsiella aerogenes test bin"));
        assertThat(bin1.getDomain(), equalTo("Bacteria"));
        assertThat(bin1.getGc(), equalTo(11));
        assertThat(bin1.getRefGenome(), equalTo("548.1022"));
        // Add it to the group.
        binGroup.addBin(bin1);
        // Create a new bin and merge the bins.
        final String node2Id = "NODE_12_length_174929_cov_82.7452_ID_1122656";
        Bin bin2 = new Bin(node2Id, 1122656, 82.7452);
        binGroup.addBin(bin2);
        binGroup.merge(bin1, bin2);
        // Verify that the bins are merged
        assertThat(bin1.isSignificant(), equalTo(true));
        assertThat(bin1.getContigs(), containsInAnyOrder(node1Id, node2Id));
        assertThat(bin1.getCoverage(), closeTo(96.107, 0.001));
        assertThat(bin1.getLen(), equalTo(2245404));
        assertThat(bin1.getTaxonID(), equalTo(548));
        assertThat(bin1.getName(), equalTo("Klebsiella aerogenes test bin"));
        assertThat(bin1.getDomain(), equalTo("Bacteria"));
        assertThat(bin1.getGc(), equalTo(11));
        assertThat(bin1.getRefGenome(), equalTo("548.1022"));
        // Verify that the contig map is updated.
        assertThat(binGroup.getContigBin(node2Id), equalTo(bin1));
        assertThat(binGroup.getContigBin(node1Id), equalTo(bin1));
        // Test JSON conversion.
        JsonObject binJson = bin1.toJson();
        Bin bin3 = new Bin(binJson);
        assertThat(bin3.isSignificant(), equalTo(true));
        assertThat(bin3.getContigs(), containsInAnyOrder(node1Id, node2Id));
        assertThat(bin3.getCoverage(), closeTo(96.107, 0.001));
        assertThat(bin3.getLen(), equalTo(2245404));
        assertThat(bin3.getTaxonID(), equalTo(548));
        assertThat(bin3.getRefGenome(), equalTo("548.1022"));
        assertThat(bin3.getName(), equalTo("Klebsiella aerogenes test bin"));
        assertThat(bin3.getDomain(), equalTo("Bacteria"));
        assertThat(bin3.getGc(), equalTo(11));
        // Test the clone check.
        assertThat(bin3.isClone(bin1), equalTo(true));
        // Try JSON conversion again with an output file.
        File outFileX = new File("data", "bin.ser");
        bin1.setOutFile(outFileX);
        bin1.close();
        binJson = bin1.toJson();
        bin3 = new Bin(binJson);
        assertThat(bin3.getOutFile().getAbsolutePath(), equalTo(outFileX.getAbsolutePath()));
        assertThat(bin3.isClone(bin1), equalTo(true));
        // Add more contigs to the group.
        File nodeFile = new File("data", "nodes.tbl");
        List<String> placedNodes = new ArrayList<String>(12);
        placedNodes.add(node1Id);
        placedNodes.add(node2Id);
        Set<String> unplacedNodes = new HashSet<String>(20);
        try (TabbedLineReader nodeStream = new TabbedLineReader(nodeFile)) {
            for (var line : nodeStream) {
                final String contigId = line.get(0);
                Bin binx = new Bin(contigId, line.getInt(1), line.getDouble(2));
                binGroup.addBin(binx);
                if (line.getFlag(3)) {
                    binGroup.merge(bin1, binx);
                    placedNodes.add(contigId);
                } else
                    unplacedNodes.add(contigId);
            }
        }
        // Insure everything that should be is in the big bin.
        var contigList = bin1.getContigs();
        assertThat(contigList.size(), equalTo(placedNodes.size()));
        for (String contigId : placedNodes)
            assertThat(contigList, hasItem(contigId));
        // There should only be one significant bin.
        var sigBins = binGroup.getSignificantBins();
        assertThat(sigBins, contains(bin1));
        // Check that the placed contigs are in the contig map correctly.
        for (String contigId : placedNodes) {
            Bin binx = binGroup.getContigBin(contigId);
            assertThat(contigId, binx, sameInstance(bin1));
        }
        // Check the unplaced nodes.
        for (String contigId : unplacedNodes) {
            Bin binx = binGroup.getContigBin(contigId);
            assertThat(binx, not(nullValue()));
            assertThat(binx.getName(), equalTo(contigId));
            assertThat(binx.isSignificant(), equalTo(false));
        }
        // Test sequence writer.
        File outDir = new File("data", "bins");
        File inFile = new File("data", "BigSample.fasta");
        binGroup.write(inFile, outDir);
        File outFile = new File(outDir, "unbinned.fasta");
        try (FastaInputStream readBackStream = new FastaInputStream(outFile)) {
            int counter = 0;
            for (Sequence seq : readBackStream) {
                String seqLabel = seq.getLabel();
                assertThat(seqLabel, unplacedNodes, hasItem(seqLabel));
                counter++;
            }
            assertThat(counter, equalTo(unplacedNodes.size()));
        }
        File binFile = new File(outDir, "bin.0001.fasta");
        try (FastaInputStream readBackStream = new FastaInputStream(binFile)) {
            int counter = 0;
            for (Sequence seq : readBackStream) {
                String seqLabel = seq.getLabel();
                assertThat(seqLabel, placedNodes, hasItem(seqLabel));
                counter++;
            }
            assertThat(counter, equalTo(placedNodes.size()));
        }
        binFile = new File(outDir, "bin.0002.fasta");
        assertThat(binFile.exists(), equalTo(false));
        // Add some counts and add the input file.
        binGroup.count("abc");
        binGroup.count("xyz", 4);
        binGroup.setInputFile(inFile);
        // Save and load the bin group.
        File saveFile = new File("data", "binGroup.ser");
        binGroup.save(saveFile);
        BinGroup loadedGroup = new BinGroup(saveFile);
        assertThat(loadedGroup.size(), equalTo(binGroup.size()));
        // Verify the contig map still works.
        for (String placed : placedNodes) {
            Bin binx = loadedGroup.getContigBin(placed);
            assertThat(placed, binx.isClone(bin1));
        }
        for (String unplaced : unplacedNodes) {
            Bin binx = loadedGroup.getContigBin(unplaced);
            Bin bin0 = binGroup.getContigBin(unplaced);
            assertThat(unplaced, binx.isClone(bin0));
        }
        // Verify the bin iteration still works.
        Map<String, Bin> originals = new HashMap<String, Bin>(binGroup.size() * 4 / 3 + 1);
        for (Bin origBin : binGroup)
            originals.put(origBin.getName(), origBin);
        for (Bin loadedBin : loadedGroup) {
            String name = loadedBin.getName();
            Bin origBin = originals.get(name);
            assertThat(name, origBin, not(nullValue()));
            assertThat(name, loadedBin.isClone(origBin));
        }
        // Verify the counts.
        assertThat(loadedGroup.getCount("abc"), equalTo(1));
        assertThat(loadedGroup.getCount("xyz"), equalTo(4));
        assertThat(loadedGroup.getCount("not-found"), equalTo(0));
        // Verify the input file.
        File loadedInFile = loadedGroup.getInputFile();
        assertThat(loadedInFile.isAbsolute(), equalTo(true));
        assertThat(loadedInFile.getAbsolutePath(), equalTo(inFile.getAbsolutePath()));
    }

    @Test
    void testMultiRefGenomes() throws IOException {
        // We need to create two bins with reference genomes, merge them, then save and load.
        // Read in some useful DNA hits.
        var hits = ProteinFinder.loadRefGenomes(new File("data", "hitList.tbl"));
        // Create a bin group.
        BinGroup binGroup = new BinGroup();
        // Create a simple bin.
        final String node1Id = "NODE_58_length_26115_cov_109.467_ID_1122748";
        Bin bin1 = new Bin(node1Id, 1122748, 109.467);
        // Make it significant.
        var hit = hits.get(node1Id);
        Genome refGenome = new Genome(new File("data", "548.1022.gto"));
        bin1.setTaxInfo(hit, "test bin", refGenome);
        // Create a new bin and merge the bins.
        final String node2Id = "NODE_12_length_174929_cov_82.7452_ID_1122656";
        Bin bin2 = new Bin(node2Id, 1122656, 82.7452);
        hit = hits.get(node2Id);
        refGenome = new Genome(new File("data", "548.1036.gto"));
        bin2.setTaxInfo(hit, "test bin 2", refGenome);
        binGroup.addBin(bin2);
        binGroup.merge(bin1, bin2);
        assertThat(bin1.getAllRefGenomes(), contains("548.1022", "548.1036"));
        // Test json conversion.
        JsonObject json = bin1.toJson();
        Bin bin3 = new Bin(json);
        assertThat(bin3.getAllRefGenomes(), contains("548.1022", "548.1036"));
        assertThat(bin3.isClone(bin1), equalTo(true));
    }



}
