package org.theseed.bins.generate;

import java.io.File;
import java.io.IOException;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.seeds.ProteinFinder;

public class FinderTest {

    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(FinderTest.class);

    /**
     * Perform a full-scale test on the finder.  This is only possible if SEEDtk_TEST_DATA has been set and
     * the "Finder" subdirectory has finder files in it.  The finder files are huge and can't be stored in GitHub.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    @Test
    public void bigFindingTest() throws IOException, InterruptedException {
        String testDirName = System.getenv("SEEDtk_TEST_DATA");
        if (! StringUtils.isBlank(testDirName)) {
            File finderDir = new File(testDirName, "Finder");
            File testFile = new File(finderDir, "PhenTrnaSyntAlph.fna");
            if (testFile.exists()) {
                log.info("Performing finder DNA species assignment test.");
                // Find the SOUR proteins in BigSample.
                var finder = new ProteinFinder(finderDir);
                File dnaFile1 = new File("data", "BigSample.fasta");
                var rolesFound = finder.findSourProteins(dnaFile1);
                // Test saving and loading the protein map.
                File protFile = new File("data", "rolesFound.ser");
                ProteinFinder.saveSeedProteins(rolesFound, protFile);
                var rolesFoundLoaded = ProteinFinder.loadSeedProteins(protFile);
                assertThat(rolesFoundLoaded.size(), equalTo(rolesFound.size()));
                for (var roleInfo : rolesFound.entrySet()) {
                    String roleId = roleInfo.getKey();
                    var locListReal = roleInfo.getValue();
                    var locListLoaded = rolesFoundLoaded.get(roleId);
                    assertThat(roleId, locListLoaded, not(nullValue()));
                    assertThat(roleId, locListLoaded.size(), equalTo(locListReal.size()));
                    for (var loc : locListLoaded) {
                        assertThat(roleId + " " + loc.toString(), locListReal, hasItem(loc));
                    }
                }
                // Now do species assignment.
                var hitsFound = finder.findRefGenomes(rolesFound, dnaFile1);
                // Test saving and loading the dna map.
                File dnaFile = new File("data", "hitsFound.ser");
                ProteinFinder.saveRefGenomes(hitsFound, dnaFile);
                var hitsFoundLoaded = ProteinFinder.loadRefGenomes(dnaFile);
                assertThat(hitsFoundLoaded.size(), equalTo(hitsFound.size()));
                for (var hitFoundInfo : hitsFound.entrySet()) {
                    String contigId = hitFoundInfo.getKey();
                    var hitLoaded = hitsFoundLoaded.get(contigId);
                    assertThat(contigId, hitLoaded, not(nullValue()));
                    assertThat(contigId, hitLoaded.isClone(hitFoundInfo.getValue()));
                }
            }

        }

    }

}
