/**
 *
 */
package org.theseed.bins.methods;

import java.io.File;
import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.bins.BinGroup;
import org.theseed.bins.BinParms;
import org.theseed.genome.Genome;
import org.theseed.sequence.seeds.ProteinFinder;

import com.github.cliftonlabs.json_simple.JsonException;

/**
 * This object represents a single phase of the binning process.  Each phase produces a terminal file that
 * indicates whether or not it is finished and provides a way to reload the binning group output by the
 * previous phase.  This allows us to restart more efficiently.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BinPhase {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinPhase.class);
    /** controlling bin processor */
    private IParms processor;

    /**
     * This interface holds the query methods we need the binning processor to suppport.
     */
    public interface IParms {

        /**
         * @return the output directory name
         */
        public File getOutDir();

        /**
         * @return the protein finder
         */
        public ProteinFinder getFinder();

        /**
         * @return the binning engine
         */
        public BinningMethod getEngine();

        /**
         * @return the binning parameters
         */
        public BinParms getParms();

        /**
         * @return a PATRIC genome with the specified ID
         *
         * @throws IOException
         */
        public Genome getGenome(String genomeId) throws IOException;

        /**
         * @return the current binning group
         */
        public BinGroup getBinGroup();

        /**
         * @return the input FASTA file for the sample
         */
        public File getInFile();

        /**
         * @return the name suffix for new genomes
         */
        public String getNameSuffix();

    }

    /**
     * Create this binning phase.
     */
    public BinPhase(IParms commandProcessor) {
        this.processor = commandProcessor;
    }

    /**
     * This method checks that the binning-group save file exists.  If it does, the phase
     * is considered complete.  This method is overridden for the reporting phase, which
     * does not modify the binning group.
     *
     * @return TRUE if this phase is already complete, else FALSE
     */
    public boolean isDone() {
        File saveFile = this.getSaveFile();
        boolean retVal = saveFile.exists();
        return retVal;
    }

    /**
     * This is a utility method used to end phases that modify the binning group.
     * The binning group is saved to the checkpoint file.
     *
     * @throws IOException
     */
    protected void finish() throws IOException {
        BinGroup binGroup = this.processor.getBinGroup();
        binGroup.save(this.getSaveFile());
    }

    /**
     * Reload the binning-group data from this phase.
     *
     * @return a binning group produced by this phase in a previous run
     *
     * @throws JsonException
     * @throws IOException
     */
    public BinGroup reload() throws IOException, JsonException {
        BinGroup retVal = new BinGroup(this.getSaveFile());
        return retVal;
    }

    /**
     * Execute this phase.
     *
     * @throws Exception
     */
    public void run() throws Exception {
        log.info("Executing {} phase.", this);
        this.runPhase();
        log.info("{} bins and unplaced contigs identified at end of {} phase.", this.processor.getBinGroup().size(), this);
        this.finish();
    }

    /**
     * Perform the tasks required to process the current phase.
     *
     * @throws Exception
     */
    protected abstract void runPhase() throws Exception;

    /**
     * Checkpoint this phase.
     *
     * @throws IOException
     */
    public void checkpoint() throws IOException {
        File outFile = new File(this.processor.getOutDir(), this.getSaveFileName());
        BinGroup binGroup = this.processor.getBinGroup();
        binGroup.save(outFile);
    }

    /**
     * @return the base name of the binning-group checkpoint file for this phase
     */
    protected abstract String getSaveFileName();

    /**
     * @return the name of the binning-group checkpoint file for this phase
     */
    protected File getSaveFile() {
        return new File(this.processor.getOutDir(), this.getSaveFileName());
    }

    @Override
    public String toString() {
        return this.getPhaseName();
    }

    /**
     * @return the name of this phase
     */
    protected abstract String getPhaseName();

    /**
     * Create a file name in the output directory.
     *
     * @param name		base name of file
     *
     * @return the desired output file name
     */
    protected File getOutFile(String name) {
        return new File(this.processor.getOutDir(), name);
    }

    /**
     * @return the master binning group
     */
    protected BinGroup getBinGroup() {
        return this.processor.getBinGroup();
    }

    /**
     * @return the input FASTA file for the sample
     */
    protected File getInFile() {
        return this.processor.getInFile();
    }

    /**
     * @return the binning parameters
     */
    protected BinParms getParms() {
        return this.processor.getParms();
    }

    /**
     * @return the protein finder
     */
    protected ProteinFinder getFinder() {
        return this.processor.getFinder();
    }

    /**
     * @return the name suffix for new genomes
     */
    protected String getNameSuffix() {
        return this.processor.getNameSuffix();
    }

    /**
     * @return the binning engine
     */
    protected BinningMethod getEngine() {
        return this.processor.getEngine();
    }

    /**
     * @return the genome with the specified ID
     *
     * @param genomeId	ID of the desired genome
     *
     * @throws IOException
     */
    protected Genome getGenome(String genomeId) throws IOException {
        return this.processor.getGenome(genomeId);
    }

    /**
     * @return the output directory name
     */
    protected File getOutDir() {
        return this.processor.getOutDir();
    }


}
