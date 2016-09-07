package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Created by rct66 on 9/6/16.
 */
@Parameters(commandDescription = "Test a model given test files.")

public class PredictionArguments {
    @Parameter(required = true, names = {"-t", "--test-set"}, description = "test set. When more than one dataset is provided (multiple -t options), " +
    "only the first will be used to produce metadata files for statistical analysis in Goby.")
    public String testSet;

    //TODO
    @Parameter(names = {"-k", "--kind"}, description = "the kind of dataset used for testing (ie test, val, training) ")
    public String type = "test";

    @Parameter(required = true, names = {"-m", "--model-path"}, description = "directory containing the model to use for prediction. ")
    public String modelPath;

    @Parameter(names = {"-l", "--model-label"}, description = "keyword specifying which version of the model to use for predictions (ie bestAUC, latest)")
    public String modelVersion = "latest";

    @Parameter(names = {"-r", "--long-report"}, description = "long report: include base count and other feature data in the prediction output")
    public boolean longReport = false;



}
