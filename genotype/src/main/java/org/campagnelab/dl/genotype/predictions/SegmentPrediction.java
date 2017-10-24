package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.performance.SegmentGenotypePredictionTest;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.util.function.Predicate;

/**
 * Represents predicted genotypes for a segment.
 */
public class SegmentPrediction extends Prediction {
    SegmentInformationRecords.ReferencePosition startPosition;
    SegmentInformationRecords.ReferencePosition endPosition;

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

    public SegmentPrediction(SegmentInformationRecords.ReferencePosition startPosition,
                             SegmentInformationRecords.ReferencePosition endPosition,
                             SegmentGenotypePrediction segmentGenotypePrediction) {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.genotypes = segmentGenotypePrediction;
    }

    public int numPredictionsWhere(SegmentGenotypePredictionTest isTrueForBaseIndex) {
        int count = 0;
        final int numBases = numBases();
        for (int baseIndex = 0; baseIndex < numBases; baseIndex++) {
            if (isTrueForBaseIndex.test(this.genotypes,baseIndex)) {
                count+=1;
            }
        }

        return count;
    }

    public int numBases() {
        return genotypes.numBases();
    }
}
