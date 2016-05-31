package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by rct66 on 5/31/16.
 */
public class QualityFeatures implements FeatureMapper {
    @Override
    public int numberOfFeatures() { return 20; }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        // do nothing. No normalization here.
    }

    int indices[] = {0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        int sumCounts = 0;
        for (BaseInformationRecords.SampleInfo sampleInfo : record.getSamplesList()) {
            for (BaseInformationRecords.CountInfo sampleCounts : sampleInfo.getCountsList()) {

                sumCounts += sampleCounts.getGenotypeCountForwardStrand();
                sumCounts += sampleCounts.getGenotypeCountReverseStrand();
            }
        }
        final int featureValue = (Math.log10(sumCounts) > featureIndex) ? 1 : 0;
        return featureValue;
    }
}
