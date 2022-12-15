/**
 *
 */
package org.theseed.bins.methods;

/**
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
    protected void runPhase() {
        // TODO code for REPORT runPhase

    }

}
