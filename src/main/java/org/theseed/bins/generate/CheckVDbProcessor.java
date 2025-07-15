/**
 *
 */
package org.theseed.bins.generate;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.ncbi.NcbiConnection;
import org.theseed.ncbi.NcbiListQuery;
import org.theseed.ncbi.NcbiTable;
import org.theseed.ncbi.XmlException;
import org.theseed.ncbi.XmlUtils;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3Connection;
import org.w3c.dom.Element;

/**
 * This command reads a brand-new CheckV database and creates the "checkv_taxon.tsv" file in the
 * "genome_db" subdirectory.  This file contains the CheckV ID, taxon ID, taxon name, and taxon string
 * for each virus type.
 *
 * The basic information is taken from the "checkv_info.tsv" file in the target directory.  This file
 * contains an ID and a taxonomy string for each CheckV virus ID.  There are two separate cases we
 * need to handle-- GenBank viruses and other viruses.  For GenBank viruses, we need to query the NCBI using
 * the accession number (which is stored in check_info as a contig ID); for the others, we need to query
 * PATRIC using the taxonomy string.
 *
 * The basic strategy is to read in the info file and collect the taxonomic strings and GenBank accession
 * strings into batches. For each taxonomy string batch, we search the P3 database for a matching taxonomic list.
 * We do this using a name lookup on the end entry. Note that because SOLR uses fuzzy matching on string fields,
 * we will need to verify the name is an exact match before using a record.  For the accession batches, we look
 * up in the NCBI Assembly table using the ASAC (assembly accession) as the key.  The records will come back in
 * a more or less random order, but we can using the Genbank Synonym to verify the accession number.  The taxonomy
 * ID is in the SpeciesTaxid field and the virus name in the Organism field.  This is a complex
 * process, which is why this command exists:  to preload everything into a flat file for later use.
 *
 * The positional parameter is the name of the checkV database directory.  This directory must contain
 * a "genome_db" subdirectory with the "checkv_info.tsv" file in it.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -b	batch size for database queries (default is 50)
 *
 *
 * @author Bruce Parrello
 *
 */
public class CheckVDbProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CheckVDbProcessor.class);
    /** connection to PATRIC */
    private P3Connection p3;
    /** map of checkV IDs to descriptors */
    private Map<String, TaxDescriptor> taxMap;
    /** map of checkV IDs to contig IDs */
    private Map<String, String> contigMap;
    /** map of taxonomic names to descriptors */
    private Map<String, TaxDescriptor> nameMap;
    /** set of saved GenBank contig IDs */
    private Set<String> genBankIds;
    /** input file name */
    private File inFileName;
    /** output file name */
    private File outFileName;
    /** connection to NCBI */
    private NcbiConnection ncbi;
    /** ncbi batch size */
    private int ncbiBatchSize;

    /**
     * This class contains the taxonomic description fields needed:  taxon ID, taxon name, and taxon string.
     * It is keyed on the lineage.
     */
    protected static class TaxDescriptor {

        /** taxon ID */
        private int taxID;
        /** taxon name */
        private String taxName;
        /** taxon lineage string */
        private String lineage;
        /** TRUE if this is a GenBank virus */
        private boolean genBankFlag;

        /**
         * Construct a taxonomic descriptor from a lineage string.
         *
         * @param lineage	source lineage string
         */
        protected TaxDescriptor(String lineage) {
            // Denote we have no taxonomic ID and this is a non-genbank virus.
            this.taxID = 0;
            this.genBankFlag = false;
            // Save the lineage.
            this.lineage = lineage;
            // Extract the bottom-level name.  If the lineage is a singleton, we use the only name.
            this.taxName = StringUtils.substringAfterLast(lineage, ";");
            if (this.taxName.isEmpty())
                this.taxName = lineage;
        }

        /**
         * Create a GenBank-style descriptor from a taxonomic ID and an organism name.
         *
         * @param taxId		species taxonomic ID
         * @param name		species organism name
         */
        protected TaxDescriptor(int taxId, String name) {
            // Initialize the taxonomic ID and denote this is a genbank virus.
            this.taxID = taxId;
            this.genBankFlag = true;
            // Strip off the parenthetical and save the name.
            this.taxName = StringUtils.removeEnd(name, " (viruses)");
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.lineage == null) ? 0 : this.lineage.hashCode());
            return result;
        }
        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof TaxDescriptor)) {
                return false;
            }
            TaxDescriptor other = (TaxDescriptor) obj;
            if (this.lineage == null) {
                if (other.lineage != null) {
                    return false;
                }
            } else if (!this.lineage.equals(other.lineage)) {
                return false;
            }
            return true;
        }

        /**
         * @return the taxonomic ID
         */
        public int getTaxID() {
            return this.taxID;
        }

        /**
         * Specify a new taxonomic ID.
         *
         * @param taxID the taxID to set
         */
        public void setTaxID(int taxID) {
            this.taxID = taxID;
        }

        /**
         * @return the taxonomic name
         */
        public String getTaxName() {
            return this.taxName;
        }

        /**
         * @return the full lineage
         */
        public String getLineage() {
            return this.lineage;
        }

        /**
         * Copy the lineage from an old descriptor.
         *
         * @param original	old descriptor
         */
        public void copyLineage(TaxDescriptor original) {
            this.lineage = original.lineage;
        }

        /**
         * @return the display name for a virus with this descriptor
         *
         * @param contigId		associated contig ID (used for non-genbank viruses)
         */
        public String getName(String contigId) {
            String retVal = this.taxName;
            if (! this.genBankFlag)
                retVal += " " + contigId;
            return retVal;
        }

    }

    // COMMAND-LINE OPTIONS

    /** number of taxon names per query batch */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "10", usage = "maximum number of taxon names to use in a query batch")
    private int batchSize;

    /** directory containing the checkV database */
    @Argument(index = 0, metaVar = "checkv_db_dir", usage = "directory containing the checkv database", required = true)
    private File dbDir;

    @Override
    protected void setDefaults() {
        this.batchSize = 50;
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        // Insure the checkV database directory exists.
        if (! this.dbDir.isDirectory())
            throw new FileNotFoundException("CheckV database directory " + this.dbDir + " does not exist or is invalid.");
        // Get the input and output file names.
        File gDbDir = new File(this.dbDir, "genome_db");
        this.inFileName = new File(gDbDir, "checkv_info.tsv");
        this.outFileName = new File(gDbDir, "checkv_taxon.tsv");
        if (! this.inFileName.canRead())
            throw new FileNotFoundException("Invalid checkV database directory:  could not find and/or read "
                    + this.inFileName + ".");
        // Verify we can write the output file.  This will fail if we can't.
        this.outFileName.createNewFile();
        // Verify the batch size.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1.");
        // Scale the NCBI batch size to be comparable to the PATRIC batch size.
        this.ncbiBatchSize = this.batchSize * 4;
        // Set up the PATRIC and NCBI connections.
        this.p3 = new P3Connection();
        this.ncbi = new NcbiConnection();
        // Create the maps.  We do an estimate based on the file size to try to optimize the size.
        int hashSize = (int) (this.inFileName.length() / 130);
        this.taxMap = new HashMap<String, TaxDescriptor>(hashSize);
        this.contigMap = new HashMap<String, String>(hashSize);
        this.nameMap = new HashMap<String, TaxDescriptor>(hashSize / 250);
        this.genBankIds = new HashSet<String>(hashSize / 2);
    }

    @Override
    protected void runCommand() throws Exception {
        // Our first task is to read the input file and create the descriptors.  For each descriptor, we will
        // add an entry to the name map record its checkv-ID and contig ID.
        log.info("Reading {}.", this.inFileName);
        int linesIn = 0;
        int badNames = 0;
        int nonViral = 0;
        try (TabbedLineReader inStream = new TabbedLineReader(this.inFileName)) {
            // Find the columns that matter.
            int idCol = inStream.findField("checkV_id");
            int contigCol = inStream.findField("contig_id");
            int lineageCol = inStream.findField("genomad_taxonomy");
            int sourceCol = inStream.findField("source");
            // Loop through the input file.
            for (var line : inStream) {
                String checkVId = line.get(idCol);
                String contigId = line.get(contigCol);
                String lineage = line.get(lineageCol);
                String source = line.get(sourceCol);
                linesIn++;
                // Connect the checkV ID to the contig ID.  The contig ID is used to form the genome name
                // in the output.
                this.contigMap.put(checkVId, contigId);
                // Is this a good lineage?
                if (! lineage.startsWith("Viruses")) {
                    // No, so we don't save the lineage and there won't be a tax ID.
                    nonViral++;
                } else {
                    // Create the taxonomic descriptor.
                    TaxDescriptor descriptor = new TaxDescriptor(lineage);
                    // Get the name.  An empty name should not happen, but might if there is a bad lineage.
                    String taxName = descriptor.getTaxName();
                    if (taxName.isBlank())
                        log.warn("Empty taxonomic name for {}.", checkVId);
                    else {
                        // Now check the name map.
                        TaxDescriptor nameDescriptor = this.nameMap.get(taxName);
                        if (nameDescriptor == null) {
                            // Here we have a new name.
                            this.nameMap.put(taxName, descriptor);
                            // Save the descriptor so we use it in our mapping later.
                            nameDescriptor = descriptor;
                        } else if (! nameDescriptor.equals(descriptor)) {
                            // Here the new name doesn't match the old one.
                            log.warn("Duplicate taxonomy for {}: \"{}\" presented but \"{}\" already found.",
                                    taxName, lineage, nameDescriptor.getLineage());
                            badNames++;
                        }
                        // If this is a genbank virus, we save the contig ID for later.
                        if (source.contentEquals("NCBI GenBank"))
                            this.genBankIds.add(contigId);
                        // Finally, connect the checkV ID to the descriptor.
                        this.taxMap.put(checkVId, nameDescriptor);
                    }
                }
            }
            log.info("Input file read:  {} lines, {} unique taxonomies, {} uncertain, {} mismatched lineages,"
                    + "{} GenBank accession IDs.", linesIn, this.nameMap.size(), nonViral, badNames,
                    this.genBankIds.size());
        }
        // Now we need to get taxon IDs for each taxonomic descriptor.
        this.processLineages();
        // The next step is to query NCBI for the taxonomic IDs of the GenBank viruses.
        this.processAccessions();
        // Now write the output.
        log.info("Writing output.");
        try (PrintWriter writer = new PrintWriter(this.outFileName)) {
            // Write the header line.
            writer.println("checkv_id\ttaxon_id\tname\ttaxonomy");
            // Loop through the contig ID map.  This will have all the checkV IDs, even the ones without
            // descriptors.
            for (var checkVEntry : this.contigMap.entrySet()) {
                String checkVId = checkVEntry.getKey();
                String contigId = checkVEntry.getValue();
                TaxDescriptor desc = this.taxMap.get(checkVId);
                // We need the taxon ID, the name, and the lineage.  These all
                // depend on whether or not there is a descriptor.
                int taxId;
                String name;
                String lineage;
                if (desc == null) {
                    taxId = 0;
                    name = contigId;
                    lineage = "";
                } else {
                    taxId = desc.getTaxID();
                    name = desc.getName(contigId);
                    lineage = desc.getLineage();
                }
                // Now print the line.
                writer.format("%s\t%d\t%s\t%s%n", checkVId, taxId, name, lineage);
            }

        }
    }

    /**
     * This method queries the NCBI to compute the species-level taxonomic ID for GenBank-sourced
     * viruses.  Our first task will be to get the assembly records from the NCBI.  This is done
     * in batches for performance.
     *
     * @throws IOException
     * @throws XmlException
     */
    private void processAccessions() throws XmlException, IOException {
        // We will create taxonomy descriptors for the NCBI contigs in here.
        var ncbiTaxMap = new HashMap<String, TaxDescriptor>(this.genBankIds.size() * 4 / 3 + 1);
        // Loop through the genbank IDs, forming batches for querying.
        log.info("{} GenBank IDs to process.", this.genBankIds.size());
        NcbiListQuery batch = new NcbiListQuery(NcbiTable.ASSEMBLY, "ASAC");
        int batchCount = 0;
        for (String contigId : this.genBankIds) {
            if (batch.size() >= this.ncbiBatchSize) {
                // The current batch is full, so query the NCBI for the tax ID and name.
                batchCount++;
                log.info("Processing NCBI batch {} with {} contig IDs.", batchCount, batch.size());
                this.processNcbiBatch(batch, ncbiTaxMap);
            }
            batch.addId(contigId);
        }
        // Process the residual batch.
        if (batch.size() > 0) {
            log.info("Processing residual NCBI batch with {} contig IDs.", batch.size());
            this.processNcbiBatch(batch, ncbiTaxMap);
        }
        // Now we run through the contig ID map updating the taxonomic descriptors for each virus for which
        // we have NCBI data.
        log.info("Updating taxonomy map with NCBI data.");
        int updateCount = 0;
        for (var contigEntry : this.contigMap.entrySet()) {
            String virusId = contigEntry.getKey();
            String contigId = contigEntry.getValue();
            TaxDescriptor newDesc = ncbiTaxMap.get(contigId);
            if (newDesc != null) {
                // Here we have NCBI info for the specified virus.  We need to store it in the main taxonomy descriptor
                // map.  First, we need the lineage from the old descriptor.
                TaxDescriptor oldDesc = this.taxMap.get(virusId);
                newDesc.copyLineage(oldDesc);
                // Store the descriptor.
                this.taxMap.put(virusId, newDesc);
                updateCount++;
            }
        }
        log.info("{} viruses updated.", updateCount);
    }

    /**
     * Process a batch of NCBI queries.  We get the summary records for the identified assemblies and
     * extract the species ID and name.
     *
     * @param batch			batch of contig IDs to process
     * @param ncbiTaxMap	output map for descriptors created
     *
     * @throws IOException
     * @throws XmlException
     */
    private void processNcbiBatch(NcbiListQuery batch, Map<String, TaxDescriptor> ncbiTaxMap) throws XmlException, IOException {
        // Run the query.
        List<Element> results = batch.run(this.ncbi);
        log.info("{} results returned from query.", results.size());
        // Now we loop through the results.  We need to extract the contig ID from the synonym list to match the records to
        // the batch.
        int badCount = 0;
        for (var result : results) {
            String contigId = null;
            Element synonyms = XmlUtils.findFirstByTagName(result, "Synonym");
            if (synonyms != null)
                contigId = XmlUtils.getXmlString(synonyms, "Genbank");
            if (contigId == null)
                badCount++;
            else {
                // Now that we have the contig ID, we can get the species information.
                int taxId = XmlUtils.getXmlInt(result, "SpeciesTaxid");
                String taxName = XmlUtils.getXmlString(result, "Organism");
                // Put it in the output hash.
                TaxDescriptor desc = new TaxDescriptor(taxId, taxName);
                ncbiTaxMap.put(contigId, desc);
            }
        }
        if (badCount > 0)
            log.warn("{} bad results from NCBI query batch.", badCount);
    }

    /**
     * This method runs through the taxonomy descriptors and uses PATRIC to fill in the taxonomic IDs
     * whenever possible.
     */
    private void processLineages() {
        // For efficiency, we break the taxon IDs into batches.
        var descBatch = new ArrayList<TaxDescriptor>(this.batchSize);
        int batchCount = 0;
        for (TaxDescriptor desc : this.nameMap.values()) {
            descBatch.add(desc);
            if (descBatch.size() >= this.batchSize) {
                batchCount++;
                log.info("Processing taxonomic-name query batch #{} with {} names.", batchCount, descBatch.size());
                this.processTaxIdBatch(descBatch);
                descBatch.clear();
            }
        }
        // Process the residual batch.
        if (descBatch.size() > 0) {
            log.info("Processing residual query batch with {} names.", descBatch.size());
            this.processTaxIdBatch(descBatch);
        }
    }

    /**
     * Query PATRIC to fill taxonomic IDs into a batch of taxonomic descriptors.
     *
     * @param descBatch		collection of descriptors to process
     */
    private void processTaxIdBatch(Collection<TaxDescriptor> descBatch) {
        // Get all the taxon names.  We use a set so we can discard records whose returned name doesn't match anything we want.
        // This happens because PATRIC does fuzzy searches for string fields.
        Set<String> taxNames = descBatch.stream().map(x -> x.getTaxName()).collect(Collectors.toSet());
        // Form a query and ask for results.
        var records = this.p3.getRecords(P3Connection.Table.TAXONOMY, "taxon_name", taxNames, "taxon_id,taxon_name");
        log.info("{} records returned from query.", records.size());
        // Now loop through the records, looking for names.
        int skipped = 0;
        int stored = 0;
        for (var record : records) {
            String taxName = KeyBuffer.getString(record, "taxon_name");
            if (! taxNames.contains(taxName))
                skipped++;
            else {
                // Here we have a real match, not a fuzzy match.  Get the descriptor.
                TaxDescriptor desc = this.nameMap.get(taxName);
                if (desc == null)
                    throw new IllegalStateException("Did not find " + taxName + " in name map.");
                else {
                    // Store the taxonomic ID.
                    desc.setTaxID(KeyBuffer.getInt(record, "taxon_id"));
                    stored++;
                }
            }
        }
        log.info("Batch processed.  {} records skipped, {} taxon IDs stored.", skipped, stored);
    }

}
