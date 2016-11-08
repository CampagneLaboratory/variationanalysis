package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV19;
import org.campagnelab.dl.varanalysis.learning.architecture.SixDenseLayersNarrower2;
import org.campagnelab.dl.varanalysis.tools.RecordingToolArguments;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Arguments for tools that train models.
 * Created by fac2003 on 8/20/16.
 */
@Parameters(commandDescription = "Train a model given training files and a validation file.")

public class TrainingArguments extends RecordingToolArguments {

    @Parameter(required = true, names = {"-t", "--training-sets"}, variableArity = true, description = "Training sets, must be provided in .sbi/.sbip format (produced with Goby3). When more than one dataset is provided (multiple -t options), the " +
            "datasets are concatenated.")
    public List<String> trainingSets = new ArrayList<>();

    @Parameter(names = {"-v", "--validation-set"}, description = "Validation set, must be provided in .parquet/.info format.")
    public String validationSet = null;


    @Parameter(names = {"-n", "--num-training"}, description = "The maximum number of training samples to train with. ")
    public int numTraining = Integer.MAX_VALUE;

    @Parameter(names = {"-x", "--num-validation"}, description = "The maximum number of validation samples to read when evaluating performance. ")
    public int numValidation = Integer.MAX_VALUE;

    @Parameter(names = {"-s", "--random-seed"}, description = "The random seed to initialize network weights. ")
    public long seed = new Date().getTime();

    @Parameter(names = "--early-stopping-num-epochs", description = "The number of epochs without performance improvement before early stopping is triggered. ")
    public int stopWhenEpochsWithoutImprovement = 10;

    @Parameter(names = {"-r", "--learning-rate"}, description = "Learning rate.")
    public double learningRate = 0.1;

    @Parameter(names = {"--dropout-rate"}, description = "Dropout rate.")
    public Double dropoutRate = null;

    @Parameter(names = "--regularization-rate", description = "Regularization rate. Disabled if not provided.")
    public Double regularizationRate =null;


    @Parameter(names = "--experimental-condition", description = "The experimental condition label used in validation loggin each epoch. ")
    public String experimentalCondition = "not_specified";

    @Parameter(names = "--mini-batch-size", description = "The size of the training minibatch")
    public int miniBatchSize = 32;

    @Parameter(names = {"--max-epochs"}, description = "The maximum number of epochs to train if early stopping does not occur")
    public int maxEpochs = Integer.MAX_VALUE;

    @Parameter(names = "--previous-model-path", description = "A model path to load parameters to continue training.")
    public String previousModelPath;

    @Parameter(names = "--previous-model-name", description = "The name of the previous model to load (i.e., \"bestAUC\", \"best\" or \"latest\" and continue training.")
    public String previousModelName = "bestAUC";

    @Parameter(names = "--net-architecture", description = "fully qualified classname that implements the choice of network architecture.")
    public java.lang.String architectureClassname = SixDenseLayersNarrower2.class.getCanonicalName();


    public String[] getTrainingSets() {
        return this.trainingSets.toArray(new String[this.trainingSets.size()]);
    }

    @Parameter(names = "--parameter-precision", description = "Parameter precision, either FP16 or FP32. Note that models trained with FP16 cannot be used on the CPU (as of DL4J 0.6.0).")
    public String precision = "FP32";

    @Parameter(names = "--feature-mapper", description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = FeatureMapperV19.class.getCanonicalName();

    @Parameter(names = {"-e", "--validate-every"}, description = "Validate only every e epochs when using early stopping. This can save time if training is much faster than evaluation.")
    public int validateEvery = 1;

    @Parameter(names = {"--error-enrichment"}, description = "When set, train with error enrichment.)")
    public boolean errorEnrichment = false;
    @Parameter(names = {"--num-errors-added"}, description = "Number of errors added to each mini-batch (only used when training with error enrichment).)")
    public int numErrorsAdded = 16;

}
