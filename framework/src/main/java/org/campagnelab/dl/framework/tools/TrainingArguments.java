package org.campagnelab.dl.framework.tools;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.campagnelab.dl.framework.tools.arguments.RecordingToolArguments;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * Arguments for tools that train models.
 * Created by fac2003 on 8/20/16.
 */
@Parameters(commandDescription = "Train a model given training files and a validation file.")

public abstract class  TrainingArguments extends RecordingToolArguments {

    @Parameter(required = true, names = {"-t", "--training-sets"}, variableArity = true, description = "Training set filenames- for example, could be provided in .sbi/.sbip format (produced with Goby3). When more than one dataset is provided (multiple -t options), the " +
            "datasets are concatenated.")
    public List<String> trainingSets = new ArrayList<>();

    @Parameter(required = true, names = {"-v", "--validation-set"}, description = "Validation set filename- for example, could be provided in .parquet/.info format.")
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
    public Double regularizationRate = null;


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
    public java.lang.String architectureClassname = defaultArchitectureClassname();

    @Parameter(names = "--ignore-cache", description = "Ignore the cache.")
    public boolean ignoreCache;

    @Parameter(names = "--memory-cache", description = "Name of the datasets to fully cache in memory. Use training,validation to " +
            "cache both training and validation set, or validation only when the training set is too large to fit fully in memory. " +
            "Can speed up training, but requires the training set to be small enough to fit in the GPU memory. The default caches" +
            "only the validation set. Use none to disable caching entirely.")
    public String memoryCache = "validation";

    @Parameter(names = "--label-smoothing-epsilon", description = "Value of epsilon for label smoothing. Zero (default) is no smoothing. Try small values (<0.1).")
    public float labelSmoothingEpsilon=0;

    public boolean memoryCacheTraining() {
        // do not cache in memory if just build a cache on disk:
        boolean result = memoryCache.length() > 1 && memoryCache.contains("training") && !buildCacheAndStop;
        if (result) {
            System.out.println("Training set will be fully cached in memory.");
        }
        return result;
    }

    public boolean memoryCacheValidation() {
        // do not cache in memory if just build a cache on disk:
        boolean result = memoryCache.length() > 1 && memoryCache.contains("validation") && !buildCacheAndStop;
        if (result) {
            System.out.println("Validation set will be fully cached in memory.");
        }
        return result;
    }

    @Parameter(names = "--add-ui-listener", description = "If true, adds a UI listener to the network to help visualize " +
            "the network and training progress. Will affect speed, so only use for tuning. ")
    public boolean addUiListener = false;

    @Parameter(names = "--ui-stats-file", description = "If not null, will save stats to this file instead of in memory" +
            " if a UI listener is being used. ")
    public String uiStatsFile = null;

    @Parameter(names = "--gpu-device", description = "Index of the GPU to use for training (0,1, up to the number of GPUs in the server).")
    public Integer deviceIndex = null;

    @Parameter(names = "--parallel", description = "When provided, trains on several GPUs in parallel.")
    public boolean parallel;

    protected abstract String defaultArchitectureClassname();

    @Parameter(names = "--build-cache-then-stop", description = "When provided, build the caches, then immediately stop.")
    public boolean buildCacheAndStop = false;

    public String[] getTrainingSets() {
        return this.trainingSets.toArray(new String[this.trainingSets.size()]);
    }

    @Parameter(names = "--parameter-precision", description = "Parameter precision, either FP16 or FP32. Note that models trained with FP16 cannot be used on the CPU (as of DL4J 0.6.0).")
    public String precision = "FP32";

    @Parameter(names = "--feature-mapper", description = "Fully qualified name of the feature mapper class.")
    public String featureMapperClassname = defaultFeatureMapperClassname();

    protected abstract String defaultFeatureMapperClassname();

    @Parameter(names = {"-e", "--validate-every"}, description = "Validate only every e epochs when using early stopping. This can save time if training is much faster than evaluation.")
    public int validateEvery = 1;

    @Parameter(names = {"--error-enrichment"}, description = "When set, train with error enrichment.)")
    public boolean errorEnrichment = false;
    @Parameter(names = {"--num-errors-added"}, description = "Number of errors added to each mini-batch (only used when training with error enrichment).)")
    public int numErrorsAdded = 16;

    @Parameter(names = {"--track"}, description = "Track either speed (SPEED) or performance (PERFS). Defaults to tracking performance metrics. Speed is useful to " +
            "optimize mini-batch-size and other factors influencing speed.")

    public TrackStyle trackingStyle = TrackStyle.PERFS;

    enum TrackStyle {
        SPEED, // show speed of each epoch in console
        PERFS // show performance metric values in console
    }

    @Parameter(names = "--eos-character", description = "If provided, use as EOS character index for alignment. If not, adds EOS to the input vocab. Must be specified if specified during pretraining, and likewise if not. ")
    public Integer eosIndex = null;
    @Parameter(names = "--previous-model-pretraining", description = "If true, previous model was pretrained, and adjust graph accordingly. ")
    public boolean previousModelPretraining = false;

    @Parameter(names = {"-amp","--advanced-model-properties"}, description = "Path to a properties file with advanced " +
            "model configuration properties. Advanced model properties are used by the computation graph architecture " +
            "in an architecture dependent manner. Property keys must follow the convention: " +
            "advancedModelProperties.classname.*, where classname is the fully qualified name of the class that needs " +
            "to access the property ")
    public File advancedModelConfiguration;

    @Parameter(names = "--reduction-rate", description = "The amount of reduction in hidden nodes to apply per layer")
    public float reductionRate = 0.36f;

    @Parameter(names = "--model-capacity", description = "A floating number that controls model capacity (i.e., number of hidden " +
            "nodes in the neural network). Use a c >=1 to control how many hidden nodes are created (#hiddenNodes=c*#inputs).")
    public float modelCapacity=1f;

    @Parameter(names = "--genomic-context-length", description = "The length of the genomic context. Must be an odd number. " +
            "When this number is smaller than the length of context included in the .sbi file, the context is trimmed while keeping the " +
            "base of interest in the center. Values between 21 and 61 are recommended." +
            " Since it is not clear what length of context provides best performance for a specific sequencing platform,  " +
            "this value should be included in a hyper-parameter search. ")
    public int genomicContextLength = Integer.MAX_VALUE;


}

