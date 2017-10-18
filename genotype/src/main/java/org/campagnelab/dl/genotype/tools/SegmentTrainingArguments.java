package org.campagnelab.dl.genotype.tools;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.mappers.SingleBaseFeatureMapperV1;

public class SegmentTrainingArguments extends TrainingArguments {
    public int numHiddenNodes=16;

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
}
