/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.theseed.bins.Bin;
import org.theseed.bins.BinGroup;

/**
 * This method produces a report on the bins.  For each bin, we output the name, the species ID,
 * the reference genome ID, the coverage, and the length.  The unplaced contigs are gathered into
 * a virtual bin to help with the reporting.  Note that the virtual bin will not be marked as
 * significant, but it provides a handy way to collect the residual contigs.
 *
 * @author Bruce Parrello
 *
 */
public class ReportPhase extends BinPhase {

    public ReportPhase(IParms commandProcessor) {
        super(commandProcessor);
    }

    @Override
    protected String getSaveFileName() {
        return "bins.json";
    }

    @Override
    protected String getPhaseName() {
        return "REPORTING";
    }

    @Override
    protected void runPhase() throws IOException {
        log.info("Gathering unplaced contigs.");
        BinGroup binGroup = this.getBinGroup();
        // The first task is to merge all the unplaced contigs into a single bin.
        Bin virtualBin = null;
        for (Bin bin : binGroup) {
            if (! bin.isSignificant()) {
                // Here we have an unplaced contig.  If it is the first, it becomes the
                // virtual bin.  Otherwise, we add it to the virtual bin.
                if (virtualBin == null) {
                    File outFile = this.getOutFile(BinGroup.UNPLACED_FASTA_NAME);
                    bin.setVirtual(outFile);
                    virtualBin = bin;
                } else
                    binGroup.merge(virtualBin, bin);
            }
        }
        // Now we produce the report on all the bins.
        File outFile = this.getOutFile("report.tbl");
        log.info("Writing summary report to {}.", outFile);
        try (PrintWriter writer = new PrintWriter(outFile)) {
            writer.println("bin_name\tspecies\tref_genome_id\tcoverage\tdna_size");
            for (Bin bin : binGroup)
                writer.format("%s\t%d\t%s\t%6.2f\t%d%n", bin.getName(), bin.getTaxonID(), bin.getRefGenome(),
                        bin.getCoverage(), bin.getLen());
        }
        // Now write out the contigs to the FASTA files.
        binGroup.write(this.getInFile(), this.getOutDir());
    }

}