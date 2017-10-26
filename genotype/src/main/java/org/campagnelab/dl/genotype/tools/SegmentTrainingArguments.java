package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.mappers.SingleBaseFeatureMapperV1;

public class SegmentTrainingArguments extends TrainingArguments {
    @Parameter(names = "--num-hidden-nodes", description = "The number of LSTM hidden nodes per layer.")
    public int numHiddenNodes=1024;
    public enum RNNKind{
        CUDNN_LSTM,
        DL4J_Graves,
        DL4J_BidirectionalGraves
    }
    @Parameter(names = "--rnn-kind", description = "Kind of RNN layer: CUDNN_LSTM, DL4J_Graves or DL4J_BidirectionalGraves.")
    public RNNKind rnnKind=RNNKind.CUDNN_LSTM;

    @Override
    protected String defaultArchitectureClassname() {
        return GenotypeSegmentsLSTM.class.getCanonicalName();
    }

    @Override
    protected String defaultFeatureMapperClassname() {

        return SingleBaseFeatureMapperV1.class.getCanonicalName();
    }

    @Parameter(names = "--num-layers", description = "The number of LSTM  layers in the model.")
    public int numLayers = 1;

    @Parameter(names = "--add-input-link", description = "Add a direct link between the input and genotype output.")
    public boolean useInputLink =false;
}
