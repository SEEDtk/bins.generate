/**
 *
 */
package org.theseed.text;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.io.template.LineTemplate;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This sub-command converts an incoming tab-delimited file to text paragraphs using a LineTemplate.
 *
 * The LineTemplate is specified using a plain text file.  The file can have multiple lines, but they
 * are concatenated into a single line.  The template contains literal text interspersed with variable
 * references and commands.  The variables are column identifiers for the main input file, while the commands
 * perform useful tasks like creating text lists or performing conditionals.  Leading and trailing whitespace
 * on each line of the template file will be trimmed.
 *
 * THe positional parameter is the name of the template file.  The input file should be on the standard input,
 * and the output will be to the standard output.
 *
 * The following command-line options are supported.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input tab-delimited file (if not the standard input)
 * -o	output text file (if not the standard output)
 *
 * @author Bruce Parrello
 *
 */
public class TemplateTextProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TemplateTextProcessor.class);
    /** compiled line template */
    private LineTemplate template;
    /** original template string */
    private String templateString;

    // COMMAND-LINE OPTIONS

    /** name of the template text file */
    @Argument(index = 0, metaVar = "templateFile.txt", usage = "name of the file containing the text of the template", required = true)
    private File templateFile;

    @Override
    protected void setPipeDefaults() {
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify that the template file exists.
        if (! this.templateFile.canRead())
            throw new FileNotFoundException("Template file " + this.templateFile + " is not found or unreadable.");
        // Read in the template.
        int count = 0;
        try (LineReader templateStream = new LineReader(this.templateFile)) {
            StringBuffer templateBuffer = new StringBuffer((int) this.templateFile.length());
            for (String templateLine : templateStream) {
                String newLine = StringUtils.strip(templateLine);
                if (templateBuffer.length() > 0 && ! newLine.startsWith("{{"))
                    templateBuffer.append(' ');
                templateBuffer.append(StringUtils.strip(templateLine));
                count++;
            }
            this.templateString = templateBuffer.toString();
        }
        log.info("{} lines read from {}.  Template string is {} characters.", count, this.templateFile, this.templateString.length());
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Compile the template.
        try {
            this.template = new LineTemplate(inputStream, this.templateString);
        } catch (ParseFailureException e) {
            throw new IOException("Parsing error in template string: " + e.getMessage());
        }
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // We basically process each input line and output the result of applying it to the template.
        int count = 0;
        long length = 0;
        long lastMessage = System.currentTimeMillis();
        log.info("Reading input file.");
        for (var line : inputStream) {
            count++;
            String translation = this.template.apply(line);
            length += translation.length();
            writer.println(translation);
            if (log.isInfoEnabled()) {
                long now = System.currentTimeMillis();
                if (now - lastMessage >= 5000)
                    log.info("{} lines read, {} characters written.", count, length);
                lastMessage = now;
            }
        }
        log.info("Using a {}-character template, {} lines were translated to {} characters of output.", this.templateString.length(),
                count, length);
    }

}
