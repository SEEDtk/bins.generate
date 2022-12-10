package org.theseed.bins.generate;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.seeds.ProteinFinder;

public class FinderTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FinderTest.class);

    @Test
    public void testFinding() throws IOException, InterruptedException {
        File testDir = new File("data", "Finder");
        var finder = new ProteinFinder(testDir);
        File dnaFile1 = new File("data", "BigSample.fasta");
        var rolesFound = finder.findSeedProteins(dnaFile1);
        assertThat(rolesFound.size(), greaterThan(0));
        var locList = rolesFound.get("PhenTrnaSyntAlph");
        assertThat(locList.size(), equalTo(2));
        List<String> nodes = locList.stream().map(x -> x.getContigId()).collect(Collectors.toList());
        assertThat(nodes, containsInAnyOrder("NODE_12_length_174929_cov_82.7452_ID_1122656",
                "NODE_36_length_51523_cov_936.232_ID_1122704"));
        File dnaFile2 = new File("data", "Magnetoc.fasta");
        rolesFound = finder.findSeedProteins(dnaFile2);
        assertThat(rolesFound.size(), greaterThan(0));
        for (var roleEntry : rolesFound.entrySet()) {
            locList = roleEntry.getValue();
           final int n = locList.size();
           log.info("{} locations found for role {}", n, roleEntry.getKey());
           int i = 0;
           while (i < n) {
               int i1 = i + 10;
               if (i1 > n) i1 = n;
               var subList = locList.subList(i, i1);
               log.info("     " + subList.stream().map(x -> x.toSproutString()).collect(Collectors.joining(", ")));
               i = i1;
           }
        }
    }

}
