/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.bins.BinGroup;
import org.theseed.bins.BinParms;
import org.theseed.bins.methods.BinningMethod;
import org.theseed.bins.methods.LoadPhase;
import org.theseed.bins.methods.ReportPhase;
import org.theseed.bins.methods.RunPhase;
import org.theseed.bins.methods.SourPhase;
import org.theseed.genome.Genome;
import org.theseed.bins.methods.BinPhase;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 *
 * This is the new binning command.  Binning involves four main phases.
 *
 * 	1.	Filtering the contigs and computing coverage to build a BinGroup.
 *  2.	Locating SOUR proteins and building the starter bins.
 *  3.	Assigning the remaining contigs to bins.
 *  4.	Writing the final report.
 *
 *  We allow multiple possible methods for the third phase.  A full binning pipeline includes assembly of reads into
 *  contigs and annotating and evaluating the bins, but this command only does the middle part, which is the binning itself.
 *
 * The positional parameters are the name of the input FASTA file, and the name of the output directory for the binning.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --clear			erase the output directory before beginning
 *
 * The following command-line options relate to the determination of the initial bins.
 *
 * --finder			name of the directory containing the protein finder files (default "Finder" in the current directory
 * --lenFilter		minimum length of a contig to be considered for the SOUR protein search (default 300)
 * --covgFilter		minimum coverage for a contig to be considered for the SOUR protein search (default 5.0)
 * --maxE			maximum e-value for BLAST hits when finding a SOUR protein (default 1e-20)
 * --refMaxE		maximum e-value for BLAST hits when finding a reference genome using a SOUR protein (default 1e-10)
 * --minLen			minimum fraction of a protein that must match in a protein BLAST hit (default 0.5)
 * --maxGap			maximum permissible gap between BLAST hits for merging (default 600)
 *
 * The following command-line options affect the entire binning process.
 *
 * --binLenFilter	minimum length of a contig to be considered for binning (default 300)
 * --binCovgFilter	minimum coverage for a contig to be considered for binning (default 5.0)
 * --xLimit			maximum run length of ambiguity characters allowed in a contig to be considered for binning
 * 					(default 10)
 * --kProt			protein kmer size (default 8)
 * --kDna			DNA kmer size (default 15)
 * --dangLen		repeat region kmer size (default 50)
 * --recipe			method to be used for binning (default STANDARD)
 * --binStrength	minimum kmer-hit differential to put a contig into a bin (default 10)
 *
 * The following command-line options relate to the PATRIC database.
 *
 * --dataAPIUrl		alternate URL for accessing the PATRIC data service used to download reference genomes
 * 					(default from environment variable "P3API_URL")
 * --nameSuffix		suffix to append to species names when computing bin names (default "clonal population")
 *
 * @author Bruce Parrello
 *
 */
public class BinProcessor extends BaseProcessor implements BinPhase.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinProcessor.class);
    /** binning parameter object */
    private BinParms parms;
    /** connection to PATRIC */
    private P3Connection p3;
    /** master bin group */
    private BinGroup binGroup;
    /** protein finder */
    private ProteinFinder finder;
    /** binning engine */
    private BinningMethod binEngine;
    /** list of phases */
    private List<BinPhase> phases;
    /** genome cache directory */
    private File genomeCacheDir;

    // COMMAND-LINE OPTIONS

    /** if specified, the output directory will be erased before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** name of the directory containing the protein-finder files */
    @Option(name = "--finder", metaVar = "FinderDir", usage = "name of the directory containing the protein-finder files")
    private File finderDir;

    /** minimum length of a contig to be considered for the SOUR protein search */
    @Option(name = "--lenFilter", metaVar = "400", usage = "minimum length of a contig to be considered for the SOUR protein search")
    private int lenFilter;

    /** minimum coverage for a contig to be considered for the SOUR protein search */
    @Option(name = "--covgFilter", metaVar = "4.0", usage = "minimum coverage for a contig to be considered for the SOUR protein search")
    private double covgFilter;

    /** maximum e-value for BLAST hits when finding a SOUR protein */
    @Option(name = "--maxE", aliases = { "-e" }, metaVar = "1e-30",
            usage = "maximum e-value for BLAST hits when finding a SOUR protein")
    private double maxE;

    /** maximum e-value for BLAST hits when finding a reference genome using a SOUR protein */
    @Option(name = "--refMaxE", metaVar = "1e-8",
            usage = "maximum e-value for BLAST hits when finding a reference genome using a SOUR protein")
    private double refMaxE;

    /** maximum permissible gap between BLAST hits for merging */
    @Option(name = "--maxGap", aliases = { "--gap", "-g" }, metaVar = "400",
            usage = "maximum permissible gap between BLAST hits for merging")
    private int maxGap;

    /** minimum fraction of a protein that must match in a protein BLAST hit */
    @Option(name = "--minLen", aliases = { "--minlen", "-l" }, metaVar = "0.80",
            usage = "minimum fraction of a protein that must match in a protein BLAST hit")
    private double minLen;

    /** minimum length of a contig to be considered for binning */
    @Option(name = "--binLenFilter", metaVar = "300", usage = "minimum length of a contig to be considered for binning")
    private int binLenFilter;

    /** minimum coverage for a contig to be considered for binning */
    @Option(name = "--binCovgFilter", metaVar = "3.0", usage = "minimum coverage for a contig to be considered for binning")
    private double binCovgFilter;

    /** maximum run length of ambiguity characters allowed in a contig to be considered for binning */
    @Option(name = "--xLimit", aliases = { "--XBad", "--scaffoldLen" }, metaVar = "50",
            usage = "maximum run length of ambiguity characters allowed in a contig to be considered for binning")
    private int xLimit;

    /** protein kmer size */
    @Option(name = "--kProt", aliases = { "--kmer", "-k" }, metaVar = "9", usage = "protein kmer size")
    private int kProt;

    /** DNA kmer size */
    @Option(name = "--kDna", metaVar = "23", usage = "DNA kmer size")
    private int kDna;

    /** repeat region kmer size */
    @Option(name = "--dangLen", metaVar = "40", aliases = { "--danglen" }, usage = "repeat region kmer size")
    private int dangLen;

    /** minimum kmer-hit differential to put a contig into a bin */
    @Option(name = "--binStrength", metaVar = "5", usage = "minimum kmer-hit differential to put a contig into a bin")
    private int binStrength;

    /** alternate URL for accessing the PATRIC data service used to download reference genomes */
    @Option(name = "--dataApiUrl", aliases = { "--dataAPIUrl" },
            usage = "alternate URL for accessing the PATRIC data service used to download reference genomes")
    private String dataApiUrl;

    /* suffix to append to species names when computing bin names */
    @Option(name = "--nameSuffix", metaVar = "\"from sample SRS100286\"",
            usage = "suffix to append to species names when computing bin names")
    private String nameSuffix;

    /** method to use for binning */
    @Option(name = "--recipe", aliases = { "--method" }, usage = "method to use for binning")
    private BinningMethod.Type methodType;

    /** name of the input FASTA file */
    @Argument(index = 0, metaVar = "sample.fasta", usage = "name of the input FASTA file")
    private File inFile;

    /** name of the output directory for the bins */
    @Argument(index = 1, metaVar = "outDir", usage = "name of the output directory for the bins")
    private File outDir;


    @Override
    protected void setDefaults() {
        this.clearFlag = false;
        this.finderDir = new File(System.getProperty("user.dir"), "Finder");
        this.binCovgFilter = 5.0;
        this.binLenFilter = 300;
        this.covgFilter = 5.0;
        this.lenFilter = 500;
        this.dangLen = 50;
        this.dataApiUrl = P3Connection.getApiUrl();
        this.kDna = 15;
        this.kProt = 8;
        this.maxE = 1e-20;
        this.refMaxE = 1e-10;
        this.maxGap = 600;
        this.minLen = 0.5;
        this.nameSuffix = "clonal population";
        this.xLimit = 30;
        this.binStrength = 10;
        this.methodType = BinningMethod.Type.STANDARD;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (! this.finderDir.isDirectory())
            throw new FileNotFoundException("Finder directory " + this.finderDir + " is not found or invalid.");
        log.info("Using finder files in {}.", this.finderDir);
        this.finder = new ProteinFinder(this.finderDir);
        // Connect to PATRIC.
        log.info("Data API URL is {}.", this.dataApiUrl);
        this.p3 = new P3Connection(this.dataApiUrl);
        // Now we process the parameters.
        this.parms = new BinParms();
        if (this.binCovgFilter < 0.0)
            throw new ParseFailureException("Binning coverage filter cannot be negative.");
        this.parms.setBinCovgFilter(this.binCovgFilter);
        if (this.binLenFilter < 0.0)
            throw new ParseFailureException("Binning length filter cannot be negative.");
        this.parms.setBinLenFilter(this.binLenFilter);
        if (this.covgFilter < 0.0)
            throw new ParseFailureException("Seed-search coverage filter cannot be negative.");
        this.parms.setCovgFilter(this.binCovgFilter);
        if (this.dangLen < 0)
            throw new ParseFailureException("Mobile-element kmer length (dangLen) cannot be negative.");
        this.parms.setDangLen(this.dangLen);
        if (this.kDna < 1)
            throw new ParseFailureException("DNA kmer length must be greater than 0.");
        this.parms.setKDna(this.kDna);
        if (this.kProt < 1)
            throw new ParseFailureException("Protein kmer length must be greater than 0.");
        this.parms.setKProt(this.kProt);
        if (this.maxE < 1e-100)
            throw new ParseFailureException("Seed-search e-value limit is too low.  Minimum is 1e-100");
        this.parms.setMaxEValue(this.maxE);
        if (this.refMaxE < 1e-100)
            throw new ParseFailureException("Reference-genome e-value limit is too low.  Minimum is 1e-100");
        this.parms.setRefMaxEValue(this.refMaxE);
        if (this.maxGap < 0)
            throw new ParseFailureException("Maximum gap size cannot be negative.");
        this.parms.setMaxGap(this.maxGap);
        if (this.minLen < 0.0 || this.minLen > 1.0)
            throw new ParseFailureException("Minimum length match fraction must be between 0 and 1.");
        this.parms.setMinLen(this.minLen);
        if (this.xLimit < 0)
            throw new ParseFailureException("Ambiguity-character limit (xLimit) cannot be negative.");
        this.parms.setXLimit(this.xLimit);
        if (this.binStrength < 1)
            throw new ParseFailureException("Bin strength cannot be less than 1.");
        this.parms.setBinStrength(this.binStrength);
        log.info("Tuning parameters are {}.", this.parms);
        // Check the name suffix.
        if (! StringUtils.isAsciiPrintable(this.nameSuffix))
            throw new ParseFailureException("Name suffix can only contain printable characters.");
        if (StringUtils.isBlank(this.nameSuffix))
            log.info("No suffix will be appended to new genome names.");
        else
            log.info("Suffix for new genome names is \"{}\".", this.nameSuffix);
        // Verify we can read the input file.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input FASTA file " + this.inFile + " is not found or unreadable.");
        // Create the binning engine.
        this.binEngine = this.methodType.create(this);
        log.info("Binning method is {}.", this.methodType);
        // Insure we have a good output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Output will be to directory {}.", this.outDir);
        // Set up the genome cache.
        this.genomeCacheDir = new File(this.outDir, "RefGenomes");
        if (! this.genomeCacheDir.isDirectory()) {
            log.info("Creating genome cache directory {}.", this.genomeCacheDir);
            FileUtils.forceMkdir(this.genomeCacheDir);
        }
        // Initialize the binning phases.
        this.phases = new ArrayList<BinPhase>();
        this.phases.add(new LoadPhase(this));
        this.phases.add(new SourPhase(this));
        this.phases.add(new RunPhase(this));
        this.phases.add(new ReportPhase(this));
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Everything is done by the phase handlers.  First, we skip the completed phases.
        int phaseIdx = 0;
        while (phaseIdx < this.phases.size() && this.phases.get(phaseIdx).isDone()) {
            log.info("{} phase is already complete for this sample.", this.phases.get(phaseIdx));
            phaseIdx++;
        }
        // Reload the last completed phase (if any).
        if (phaseIdx == 0)
            this.binGroup = new BinGroup();
        else {
            BinPhase previous = this.phases.get(phaseIdx - 1);
            log.info("Reloading binning data from {} phase.", previous);
            this.binGroup = previous.reload();
        }
        // Run the remaining phases.
        while (phaseIdx < this.phases.size()) {
            this.phases.get(phaseIdx).run();
            phaseIdx++;
        }
    }

    @Override
    public File getOutDir() {
        return this.outDir;
    }

    @Override
    public ProteinFinder getFinder() {
        return this.finder;
    }

    @Override
    public BinningMethod getEngine() {
        return this.binEngine;
    }

    @Override
    public BinParms getParms() {
        return this.parms;
    }

    @Override
    public Genome getGenome(String genomeId) throws IOException {
        // Different binning methods require different parts of the genome, so we just get everything.
        return P3Genome.load(p3, genomeId, Details.FULL, genomeCacheDir);
    }

    @Override
    public BinGroup getBinGroup() {
        return this.binGroup;
    }

    @Override
    public File getInFile() {
        return this.inFile;
    }

    @Override
    public String getNameSuffix() {
        return this.nameSuffix;
    }

}
