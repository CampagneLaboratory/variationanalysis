package org.campagnelab.dl.somatic.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for the Show tool.
 */

@Parameters(commandDescription = "Show records of a dataset. The dataset must be provided in the .sbi/sbip format.")

public class ShowArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--dataset"}, description = "Path to the dataset to show.")
    public String datasetFilename;

    @Parameter(required = false, names = {"-r", "--report"}, description = "Type of report, one of POSITIONS, PROTOBUFF")
    ShowReportTypes reportType;

    @Parameter(names = {"-n", "--num-examples"}, description = "Number of examples to show.")
    public int showN = Integer.MAX_VALUE;

    @Parameter(names = {"-p", "--predictions"}, description = "BinaryClassPrediction file with indices to show. Set to - to read from standard input.")
    public String predictionFilter = null;

    public enum ShowReportTypes {
        PROTOBUFF,
        POSITIONS
    }
}
