/**
 *
 */
package org.theseed.bins;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;
import org.theseed.counters.CountMap;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * @author Bruce Parrello
 *
 */
class TestContigFilter {


    @Test
    void testContigFiltering() throws IOException {
        File testFasta = new File("data", "testSample.fasta");
        // We use the default parms with a tweak to bin-coverage so we can get a normal-low-coverage case.
        BinParms parms = new BinParms().setBinCovgFilter(4.0);
        ContigFilter filter = new ContigFilter(parms);
        // We track counters in here.
        CountMap<String> counters = new CountMap<String>();
        try (FastaInputStream inStream = new FastaInputStream(testFasta)) {
            // We run each sequence individually, insuring the results are what we expect.
            Sequence seq = inStream.next();
            Bin bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.BAD));
            assertThat(counters.getCount("contig-bad-AmbiguityRun"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(5.4, 0.1));
            seq = inStream.next();
            bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.BAD));
            assertThat(counters.getCount("contig-bad-TooShort"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(6, 0.1));
            seq = inStream.next();
            bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.SEED_USABLE));
            assertThat(counters.getCount("contig-seed-usable"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(106, 0.1));
            seq = inStream.next();
            bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.NORMAL));
            assertThat(counters.getCount("contig-normal-TooShort"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(10.5, 0.1));
            seq = inStream.next();
            bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.BAD));
            assertThat(counters.getCount("contig-bad-LowCoverage"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(1.065, 0.001));
            seq = inStream.next();
            bin = filter.computeBin(seq, counters);
            assertThat(bin.getStatus(), equalTo(Bin.Status.NORMAL));
            assertThat(counters.getCount("contig-bad-LowCoverage"), equalTo(1));
            assertThat(bin.getContigs(), contains(seq.getLabel()));
            assertThat(bin.getLen(), equalTo(seq.length()));
            assertThat(bin.getCoverage(), closeTo(4.778013, 0.00001));
            // Now read the rest and check the counters at the end.
            while (inStream.hasNext()) {
                seq = inStream.next();
                bin = filter.computeBin(seq, counters);
                assertThat(bin.getContigs(), contains(seq.getLabel()));
                assertThat(seq.getLabel(), bin.getLen(), equalTo(seq.length()));
            }
            assertThat(counters.getCount("contig-bad-AmbiguityRun"), equalTo(1));
            assertThat(counters.getCount("contig-bad-TooShort"), equalTo(3));
            assertThat(counters.getCount("contig-bad-LowCoverage"), equalTo(5));
            assertThat(counters.getCount("contig-normal-TooShort"), equalTo(4));
            assertThat(counters.getCount("contig-seed-usable"), equalTo(2));
            assertThat(counters.getCount("contig-normal-LowCoverage"), equalTo(6));
        }
    }

    @Test
    void testBinGroupLoad() throws IOException {
        final Pattern NORMAL_PATTERN = Pattern.compile("NODE_[346]\\d.*");
        // Set up to do a load.  We expect two records in the reduce file (Nodes 30 and 31) and
        // 12 in the group (Nodes 60 through 65, 40 through 43, and 30 and 31).
        File inFile = new File("data", "testSample.fasta");
        File reduceFile = new File("data", "reduced.ser");
        BinParms parms = new BinParms().setBinCovgFilter(4.0);
        BinGroup group = new BinGroup(inFile, parms, reduceFile);
        // Verify the reduced file.
        try (FastaInputStream inStream = new FastaInputStream(reduceFile)) {
            Sequence seq = inStream.next();
            assertThat(seq.getLabel(), equalTo("NODE_30_coverage_106"));
            seq = inStream.next();
            assertThat(seq.getLabel(), equalTo("NODE_31_coverage_106"));
            assertThat(inStream.hasNext(), equalTo(false));
        }
        // Verify the bin group.
        assertThat(group.size(), equalTo(12));
        for (Bin bin : group) {
            var contigs = bin.getContigs();
            String name = bin.getName();
            assertThat(name, bin.getStatus(), not(equalTo(Bin.Status.BAD)));
            assertThat(name, contigs.size(), equalTo(1));
            assertThat(name, NORMAL_PATTERN.matcher(contigs.get(0)).matches());
        }
    }

}
