package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for the Show tool.
 */

@Parameters(commandDescription = "Show records of a dataset. The dataset must be provided in the .sbi/sbip format.")

public abstract class ShowArguments  implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--dataset"}, description = "Path to the dataset to show.")
    public String datasetFilename;

    @Parameter(required = false, names = {"-r", "--report"}, description = "Type of report, one of POSITIONS, PROTOBUFF, COUNTS, FREQ_COUNTS")
    public String  reportType=defaultReportType();

    protected abstract String defaultReportType();

    @Parameter(names = {"-n", "--num-examples"}, description = "Number of examples to show.")
    public int showN = Integer.MAX_VALUE;

    @Parameter(required = true, names = {"-m", "--model-path"}, description = "directory containing the model to use for prediction. " +
            "The prediction output will be stored in this directory following the pattern <model-label>-<kind>.tsv")
    public String modelPath;

    @Parameter(names = {"-l", "--model-name"}, description = "keyword specifying which specific model to use for predictions (ie bestAUC, latest)")
    public String modelName = "bestAUC";

    @Parameter(names = {"-p", "--predictions"}, description = "BinaryClassPrediction file with indices to show. Set to - to read from standard input.")
    public String predictionFilter = null;



}
