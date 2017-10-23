package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * Represents predicted genotypes for a segment.
 */
public class SegmentPrediction extends Prediction {
    String chromosome;
    int position;

    SegmentGenotypePrediction genotypes;

    // predicted colors, one per allele ( up to ploidy), predictedColors[alleleIndex][baseIndex].
    // note that alleleIndex identifies the allele in the predicted genotype. 0 is first allele in 0/1/2,
    // 1 second, and so on.
    //byte[][] predictedColors;

    /**
     * Meta data, either populated from the segment record, or from meta-data in the cache.
     */
    // isVariant[baseIndex] is true when the base does not match the reference.
    boolean[] isVariant;
    // isIndel[baseIndex] is true when the base participates to an indel.
    boolean[] isIndel;
}
