package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.mappers.SingleBaseFeatureMapperV1;

public class SegmentTrainingArguments extends TrainingArguments {
    @Override
    protected String defaultArchitectureClassname() {
        return GenotypeSegmentsLSTM.class.getCanonicalName();
    }

    @Override
    protected String defaultFeatureMapperClassname() {

        return SingleBaseFeatureMapperV1.class.getCanonicalName();
    }
}
