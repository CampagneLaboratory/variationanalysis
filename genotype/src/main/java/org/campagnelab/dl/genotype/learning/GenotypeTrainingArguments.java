package org.campagnelab.dl.genotype.learning;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV1;

/**
 * Arguments specific to somatic model training.
 */
public class GenotypeTrainingArguments extends TrainingArguments {

    @Parameter(names = "--auc-clip-max-observations", description = "The maximum number of observations to sample when evaluating the AUC. ")
    public int aucClipMaxObservations = 10000;

    @Override
    protected String defaultArchitectureClassname() {
        return GenotypeSixDenseLayersNarrower2.class.getCanonicalName();
    }

    @Override
    protected String defaultFeatureMapperClassname() {
        return GenotypeMapperV1.class.getCanonicalName();
    }
}
