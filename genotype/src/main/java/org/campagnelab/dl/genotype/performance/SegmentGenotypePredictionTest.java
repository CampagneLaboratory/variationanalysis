package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.genotype.predictions.SegmentGenotypePrediction;

public interface SegmentGenotypePredictionTest {
    boolean test(SegmentGenotypePrediction segmentGenotypePrediction, Integer baseIndex);
}
