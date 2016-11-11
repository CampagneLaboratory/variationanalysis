package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Features that encode the magnitude of a number. A feature is 1 when the log10 of the number is larger than the feature index.
 * 31 features are sufficient to encode the magnitude of the largest integer.
 * <p>
 * Created by fac2003 on 5/24/16.
 */
public class MagnitudeFeatures extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>{
    @Override
    public int numberOfFeatures() {
        return 31;
    }

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
