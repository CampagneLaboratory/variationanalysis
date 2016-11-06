package org.campagnelab.dl.varanalysis.learning;

import com.beust.jcommander.Parameter;

/**
 * Arguments specific to somatic model training.
 */
public class SomaticTrainingArguments extends TrainingArguments {
    @Parameter(names = "--trio", description = "Use to train trio models. The training and validation datasets must have three samples, parents first, patient last.")
    public boolean isTrio = false;
    @Parameter(names = "--auc-clip-max-observations", description = "The maximum number of observations to sample when evaluating the AUC. ")
    public int aucClipMaxObservations = 10000;

}
