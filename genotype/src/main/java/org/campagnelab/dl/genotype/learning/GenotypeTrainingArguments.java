package org.campagnelab.dl.genotype.learning;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV1;

/**
 * Arguments specific to somatic model training.
 */
public class GenotypeTrainingArguments extends TrainingArguments {

    @Parameter(names = "--early-stopping-measure", description = "Name of the measure to monitor to stop early stopping. " +
            "One of score, AUC+F1, F1, Concordance, Accuracy, AUC (of correct variant identification), " +
            "Recall, Precision.")
    public String earlyStoppingMeasureName="score";

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

    @Parameter(names = "--ploidy", description = "The organism ploidy (2 for humans, more for some plants). ")
    public int ploidy=2;
}
