package org.campagnelab.dl.somatic.learning;

import com.beust.jcommander.Parameter;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2WithFrequencyAndBase;
import org.campagnelab.dl.somatic.mappers.FeatureMapperV25;

/**
 * Arguments specific to somatic model training.
 */
public class SomaticTrainingArguments extends TrainingArguments {
    @Parameter(names = "--trio", description = "Use to train trio models. The training and validation datasets must have three samples, parents first, patient last.")
    public boolean isTrio = false;
    @Parameter(names = "--auc-clip-max-observations", description = "The maximum number of observations to sample when evaluating the AUC. ")
    public int aucClipMaxObservations = 10000;

    @Parameter(names = "--early-stopping-measure", description = "Name of the measure to monitor to stop early stopping. One of score or AUC.")
    public String earlyStoppingMeasureName="AUC";

    @Override
    protected String defaultArchitectureClassname() {
        return SixDenseLayersNarrower2WithFrequencyAndBase.class.getCanonicalName();
    }

    @Override
    protected String defaultFeatureMapperClassname() {
        return FeatureMapperV25.class.getCanonicalName();
    }
}
