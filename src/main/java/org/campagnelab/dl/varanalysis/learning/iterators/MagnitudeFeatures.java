package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.format.PosRecord;
import org.campagnelab.dl.varanalysis.format.SampleRecord;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Features that encode the magnitude of a number. A feature is 1 when the log10 of the number is larger than the feature index.
 * 31 features are sufficient to encode the magnitude of the largest integer.
 *
 * Created by fac2003 on 5/24/16.
 */
public class MagnitudeFeatures implements FeatureMapper {
    @Override
    public int numberOfFeatures() {
        return 31;
    }

    int indices[] = {0, 0};

    @Override
    public void mapFeatures(PosRecord record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    @Override
    public float produceFeature(PosRecord record, int featureIndex) {
        int sumCounts = 0;
        for (SampleRecord sample : record.getSamples()) {
            for (int value : sample.getCounts()) {
                sumCounts += value;
            }
        }
        return (Math.log10(sumCounts) > featureIndex) ? 1 : 0;
    }
}
