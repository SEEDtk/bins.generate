/**
 *
 */
package org.theseed.text;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.FieldInputStream;
import org.theseed.io.LineReader;
import org.theseed.io.template.LineTemplate;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This sub-command converts incoming files to text paragraphs using a LineTemplate.
 *
 * The LineTemplate is specified using a plain text file.  We allow multiple templates, some of which are primary
 * and some of which can be linked.
 *
 * Each template can have multiple lines, but they are concatenated into a single line.  The template contains
 * literal text interspersed with variable references and commands.  The variables are column identifiers
 * for the relevant input file, while the commands perform useful tasks like creating text lists or performing
 * conditionals.  Leading and trailing whitespace on each line of the template file will be trimmed.  If the
 * first non-white character in a line is part of a literal, then a space is added in front when it is concatenated.
 *
 * The template file generally contain multiple templates for different files.  Each template corresponds to a file
 * name in the positional parameters, in order.  A main template has a header that simply says "#main".  These
 * templates generate output.  There are also linked templates. These begin with a header line that says "#linked",
 * and represent data lines that are joined with a main file using a column in that file and a column in the linked
 * file.  The column specifiers are specified as positional parameters on the header line, space-delimited. So,
 *
 * 		#linked patric_id feature_id
 *
 * would join based on the "patric_id" field in the main file using the "feature_id" field in the secondary file.  For
 * each matching line in the secondary file, the applied template text will be added to the current-line output of
 * the main file.  The linked files must follow the main file to which they apply.
 *
 * The first record of the template file absolutely must be "#main".  This is how we know we have the right kind of
 * file.  Subsequent lines that begin with "##" are treated as comments.
 *
 * THe positional parameters are the name of the template file and the name of the input files that correspond one-to-one
 * with the templates.  The output will be to the standard output.
 *
 * The input files are all field-input streams and the type is determined by the filename extension.  An extension of
 * ".tbl", ".tab", ".txt", or ".tsv" implies a tab-delimited file and an extension of ".json" a JSON list file.
 *
 * The following command-line options are supported.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output text file (if not the standard output)
 *
 * @author Bruce Parrello
 *
 */
public class TemplateTextProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TemplateTextProcessor.class);
    /** compiled main line template */
    private LineTemplate template;
    /** list of linked-template descriptors */
    private List<LinkedTemplateDescriptor> linkedTemplates;
    /** index number of current template, to match to input file index */
    private int fileIdx;

    // COMMAND-LINE OPTIONS

    /** name of the template text file */
    @Argument(index = 0, metaVar = "templateFile.txt", usage = "name of the file containing the text of the template", required = true)
    private File templateFile;

    /** name of the primary input file */
    @Argument(index = 1, metaVar = "inputFile1.tbl inputFile2.tbl", usage = "names of the data input files", required = true)
    private List<File> inputFiles;

    @Override
    protected void setReporterDefaults() {
        this.inputFiles = new ArrayList<File>();
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that the template file exists.
        if (! this.templateFile.canRead())
            throw new FileNotFoundException("Template file " + this.templateFile + " is not found or unreadable.");
        // Insure that all the input files exist.
        for (File inputFile : inputFiles) {
            if (! inputFile.canRead())
                throw new FileNotFoundException("Input file " + inputFile + " is not found or unreadable.");
        }
        // Initialize the linking structures.
        this.linkedTemplates = new ArrayList<LinkedTemplateDescriptor>(this.inputFiles.size());
        this.fileIdx = 0;
        this.template = null;
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // We start by reading a main template, then all its linked templates.  When we hit end-of-file or a
        // #main marker, we run the main template and output the results.
        try (LineReader templateStream = new LineReader(this.templateFile)) {
            // We will buffer each template group in here.  A group starts with a #main header and runs through
            // the next #main or end-of-file.
            List<String> templateGroup = new ArrayList<String>(100);
            // Special handling is required for the first header.
            Iterator<String> streamIter = templateStream.iterator();
            if (! streamIter.hasNext())
                throw new IOException("No data found in template file.");
            String savedHeader = streamIter.next();
            if (! StringUtils.startsWith(savedHeader, "#main"))
                throw new IOException("Template file does not start with #main header.");
            // Now we have the main header saved, and we can process each template group.
            // Loop through the template lines.
            for (var templateLine : templateStream) {
                if (StringUtils.startsWith(templateLine, "#main")) {
                    // New group starting.  Process the old group.
                    this.processGroup(savedHeader, templateGroup, writer);
                    // Set up for the next group.
                    templateGroup.clear();
                    savedHeader = templateLine;
                } else if (! templateLine.startsWith("##"))
                    templateGroup.add(templateLine);
            }
            // Process the residual group.
            this.processGroup(savedHeader, templateGroup, writer);
        }
    }

    /**
     * Process a single template group.  The first template is the main one, and the
     * remaining templates are all linked.
     *
     * @param savedHeader		saved header line for the main template
     * @param templateGroup		list of non-comment lines making up the template group
     * @param writer			output print writer
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    private void processGroup(String savedHeader, List<String> templateGroup, PrintWriter writer)
            throws IOException, ParseFailureException {
        if (templateGroup.isEmpty())
            log.warn("Empty template group skipped.");
        else {
            // We will buffer each template's data lines in here.
            List<String> templateLines = new ArrayList<String>(templateGroup.size());
            Iterator<String> groupIter = templateLines.iterator();
            // Process the main template.
            String linkHeader = null;
            while (groupIter.hasNext() && linkHeader == null) {
                String templateLine = groupIter.next();
                if (templateLine.startsWith("#linked"))
                    linkHeader = templateLine;
                else
                    templateLines.add(templateLine);
            }
            if (templateLines.isEmpty())
                throw new IOException("Template group has no main template.");
            // Get the input file for the main template.
            File mainFile = this.getSourceFile();
            try (FieldInputStream mainStream = FieldInputStream.create(mainFile)) {
                this.buildMainTemplate(mainStream, savedHeader, templateLines);
                // If there are any linked templates, we build them here.
                if (linkHeader != null) {
                    templateLines.clear();
                    while (groupIter.hasNext()) {
                        String templateLine = groupIter.next();
                        if (templateLine.startsWith("#linked")) {
                            // We have a new template, so build the one we've accumulated.
                            this.buildLinkedTemplate(mainStream, linkHeader, templateLines);
                            // Set up for the next template.
                            templateLines.clear();
                            linkHeader = templateLine;
                        } else
                            templateLines.add(templateLine);
                    }
                    // Build the residual template.
                       this.buildLinkedTemplate(mainStream, linkHeader, templateLines);
                    // Now link up the keys to the main file.
                    log.info("Finding keys for {} linked templates.", this.linkedTemplates.size());
                    for (var linkedTemplate : this.linkedTemplates)
                        linkedTemplate.findMainKey(mainStream);
                }
                // All the templates are compiled.  Now we run through the main file.
                // We want to count the number of lines read, the number of linked lines added, and the total text length generated.
                int count = 0;
                int linked = 0;
                long length = 0;
                long lastMessage = System.currentTimeMillis();
                // This list is used to buffer the main template and the linked ones.
                List<String> translations = new ArrayList<String>(this.linkedTemplates.size() + 1);
                log.info("Reading input file.");
                for (var line : mainStream) {
                    count++;
                    String translation = this.template.apply(line);
                    // Only produce output if the result is nonblank.
                    if (! StringUtils.isBlank(translation)) {
                        translations.add(translation);
                        // Get the linked templates.
                        for (var linkedTemplate : this.linkedTemplates) {
                            List<String> found = linkedTemplate.getStrings(line);
                            linked += found.size();
                            translations.addAll(found);
                        }
                        // Join them all together.
                        translation = StringUtils.join(translations, ' ');
                        // Now print the result.
                        length += translation.length();
                        writer.println(translation);
                        if (log.isInfoEnabled()) {
                            long now = System.currentTimeMillis();
                            if (now - lastMessage >= 5000)
                                log.info("{} lines read, {} characters written.", count, length);
                            lastMessage = now;
                        }
                        translations.clear();
                    }
                }
                if (this.linkedTemplates.size() > 0)
                    log.info("{} linked lines were incorporated from {} templates.", linked, this.linkedTemplates.size());
                log.info("{} lines were translated to {} characters of output.", count, length);
            }
        }
    }

    /**
     * Build a main template.  The main template clears the linked-template queue and is stored in the
     * main data structures.
     *
     * @param mainStream		input file stream to which the template will be applied
     * @param savedHeader		header line for this linked template
     * @param templateLines		list of template text lines
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private void buildMainTemplate(FieldInputStream mainStream, String savedHeader, List<String> templateLines)
            throws IOException, ParseFailureException {
        // Form the template into a string.
        String templateString = LinkedTemplateDescriptor.buildTemplate(templateLines);
        // Compile the string.
        log.info("Compiling main template from {} lines.", templateLines.size());
        this.template = new LineTemplate(mainStream, templateString);
        // Denote there are no linked templates yet.
        this.linkedTemplates.clear();
    }

    /**
     * Build a linked template.  The saved header must be parsed to extract the two key field names.
     * Then the template descriptor is built from the key field names, the template lines, and the
     * linked file name.  Finally, the linked template is added to the linked-template queue.
     *
     * @param mainStream		input file stream for the main template
     * @param savedHeader		header line for this linked template
     * @param templateLines		list of template text lines
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private void buildLinkedTemplate(FieldInputStream mainStream, String savedHeader, List<String> templateLines)
            throws ParseFailureException, IOException {
        // Insure we have data in the template.
        if (templateLines.isEmpty())
            throw new IOException("Empty linked template with header \"" + savedHeader + "\".");
        // Parse the header.
        String[] tokens = StringUtils.split(savedHeader);
        String mainKey;
        String linkKey;
        switch (tokens.length) {
        case 1 :
            throw new ParseFailureException("Template header \"" + savedHeader + "\" has two few parameters.");
        case 2 :
            // Single key name, so it is the same for both files.
            mainKey = tokens[1];
            linkKey = tokens[1];
            break;
        default :
            // Two key names, so we use both.
            mainKey = tokens[1];
            linkKey = tokens[2];
            break;
        }
        File linkFile = this.getSourceFile();
        // Create the template and add it to the queue.
        var template = new LinkedTemplateDescriptor(mainKey, linkKey, templateLines, linkFile);
        this.linkedTemplates.add(template);
    }

    /**
     * @return the next available source file, to use for the current template
     *
     * @throws ParseFailureException
     */
    private File getSourceFile() throws ParseFailureException {
        // Get the template input source file.
        if (this.fileIdx >= this.inputFiles.size())
            throw new ParseFailureException("Too input files in parameter list for number of templates specified.");
        File retVal = this.inputFiles.get(this.fileIdx);
        this.fileIdx++;
        return retVal;
    }

}
