package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.ToolArguments;

/**
 * Arguments for the Predict Tool.
 */
@Parameters(commandDescription = "Perform model prediction on a dataset. Dataset must be provided in the .sbi/sbip format.")

public class PredictArguments implements ToolArguments {
    @Parameter(required = true, names = {"-i", "--dataset"}, description = "Path to the dataset to where prediction will be evaluated.")
    public String testSet;

    @Parameter(names = {"-k", "--kind"}, description = "Kind of dataset used for testing (ie test, validation, training). The kind of dataset is used to construct the output filename. ")
    public String type = "test";

    @Parameter(required = true, names = {"-m", "--model-path"}, description = "directory containing the model to use for prediction. " +
            "The prediction output will be stored in this directory following the pattern <model-label>-<kind>.tsv")
    public String modelPath;

    @Parameter(names = {"-l", "--model-name"}, description = "keyword specifying which specific model to use for predictions (ie bestAUC, latest)")
    public String modelName = "latest";

    @Parameter(names = {"-r", "--long-report"}, description = "long report: include base count and other feature data in the prediction output")
    public boolean longReport = false;

    @Parameter(names = {"-n", "--num-examples"}, description = "number of examples to iterate over in the test set. useful for quickly approximating performance with fewer examples.")
    public int scoreN = Integer.MAX_VALUE;

    @Parameter(names = {"--minibatch-size"}, description = "Number of records in minibatch.")
    public int miniBatchSize = 512;

    @Parameter(names = {"--records-for-auc"}, description = "Number of records to use when evaluating AUC. Precision increases with larger values, but calculation is O(n^2) on this number.")
    public int numRecordsForAUC = 50000;

    @Parameter(names = {"--correctness-filter"},
            description = "When provided, filter output by correctness. For instance --correctness-filter wrong will only print wrong predictions. Alternatively --correctness-filter correct prings only correct predictions. ")
    public String correctnessFilter = null;
    @Parameter(names = {"--filter-p-min"},
            description = "Only output prediction with Math.max(pLabelTrue,pLabelFalse)>= x, --filter-p-min x.")
    public double pFilterMinimum = 0;
    @Parameter(names = {"--filter-p-max"},
            description = "Only output prediction with Math.max(pLabelTrue,pLabelFalse)<= x, --filter-p-max x.")
    public double pFilterMaximum = 1;

    @Parameter(names = {"-f", "--to-file"}, description = "Write output to a file. If not provided, write to stdout.")
    public boolean toFile = false;

    @Parameter(names = { "--filter-auc-observations"}, description = "When true, estimate AUC only from the filtered observations. ")
    public boolean filterAucObservations;

}
