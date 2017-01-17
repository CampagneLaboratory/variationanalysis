package org.campagnelab.dl.genotype.learning;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV1;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV24;

/**
 * Arguments specific to somatic model training.
 */
public class GenotypeTrainingArguments extends TrainingArguments {

    @Parameter(names = "--early-stopping-measure", description = "Name of the measure to monitor to stop early stopping. " +
            "One of score, AUC_VxR, F1, Concordance, Accuracy, AUC (of correct variant identification), " +
            "Recall, Precision.")
    public String earlyStoppingMeasureName = "score";

    @Parameter(names = "--auc-clip-max-observations", description = "The maximum number of observations to sample when evaluating the AUC. ")
    public int aucClipMaxObservations = 10000;

    @Parameter(names = "--variant-loss-weight", description = "The weight of variants in the loss function. " +
            "A number larger than 1 gives more weight to variant sites than non-variant sites and helps drive optimization" +
            "towards models that do well for variant prediction. This value should be included in a hyper-parameter search. ")
    public double variantLossWeight = 50;
    @Parameter(names = "--genomic-context-length", description = "The length of the genomic context. Must be an odd number. " +
            "When this number is smaller than the length of context included in the .sbi file, the context is trimmed while keeping the " +
            "base of interest in the center. Values between 21 and 61 are recommended." +
            " Since it is not clear what length of context provides best performance for a specific sequencing platform,  " +
            "this value should be included in a hyper-parameter search. ")
    public int genomicContextLength = Integer.MAX_VALUE;

    @Parameter(names = "--lstm-feature-mapper", description = "Fully qualified classname of feature mapper for LSTM used to train indels")
    public String lstmFeatureMapperClassname = null;

    @Parameter(names = "--num-layers", description = "The number of dense layers in the feedforward model.")
    public int numLayers = 5;

    @Parameter(names = "--num-lstm-layers",
            description = "The number of LSTM hidden layers, if an LSTM is being used for indels. ")
    public int numLSTMLayers = 3;

    @Parameter(names = "--num-pre-vertex-layers", description = "The number of dense layers in the feedforward model before the merge vertex. ")
    public int numPreVertexLayers = 3;

    @Parameter(names = "--num-lstm-outputs", description = "The number of outputs from the last layer in the LSTM")
    public int numLSTMOutputs = 18;

    @Parameter(names = "--model-capacity", description = "A floating number that controls model capacity (i.e., number of hidden " +
            "nodes in the neural network). Use a c >=1 to control how many hidden nodes are created (#hiddenNodes=c*#inputs).")
    public float modelCapacity=1f;

    @Override
    protected String defaultArchitectureClassname() {
        return GenotypeSixDenseLayersNarrower2.class.getCanonicalName();
    }

    @Override
    protected String defaultFeatureMapperClassname() {
        return GenotypeMapperV24.class.getCanonicalName();
    }

    @Parameter(names = "--ploidy", description = "The organism ploidy (2 for humans, more for some plants). ")
    public int ploidy = 2;
}
