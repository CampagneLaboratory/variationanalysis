package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * Represent the "genotype" output predictions in human readable form.
 *
 */
class SegmentGenotypePrediction extends Prediction {
    // predictedGenotypes[siteIndex] contains the predicted genotype (e.g., G/-, T/T, etc.)
    String[] predictedGenotypes;

    // trueGenotypes[siteIndex] contains the true genotype (e.g., G/-, T/T, etc.), obtained from the segment record.
    String[] trueGenotypes;

    // probability of prediction, one per base predicted
    float[] probabilities;
}