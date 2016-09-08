package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Created by rct66 on 9/6/16.
 */
@Parameters(commandDescription = "Predict the label in a dataset. Dataset must be provided in the .sbi/sbip format.")

public class PredictionArguments {
    @Parameter(required = true, names = {"-i", "--dataset"}, description = "Path to the dataset to where prediction will be evaluated.")
    public String testSet;

    //TODO
    @Parameter(names = {"-k", "--kind"}, description = "Kind of dataset used for testing (ie test, validation, training). The kind of dataset is used to construct the output filename. ")
    public String type = "test";

    @Parameter(required = true, names = {"-m", "--model-path"}, description = "directory containing the model to use for prediction. " +
            "The prediction output will be stored in this directory following the pattern <model-label>-<kind>.tsv")
    public String modelPath;

    @Parameter(names = {"-l", "--model-label"}, description = "keyword specifying which version of the model to use for predictions (ie bestAUC, latest)")
    public String modelVersion = "latest";

    @Parameter(names = {"-r", "--long-report"}, description = "long report: include base count and other feature data in the prediction output")
    public boolean longReport = false;

    @Parameter(names = {"-n", "--num-examples"}, description = "number of examples to iterate over in the test set. useful for quickly approximating performance with fewer examples.")
    public int scoreN = Integer.MAX_VALUE;




}
