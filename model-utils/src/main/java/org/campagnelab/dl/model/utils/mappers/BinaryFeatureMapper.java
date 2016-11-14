package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Maps an integer to 32 binary features, using one hot encoding.
 * Created by fac2003 on 7/12/16.
 */
public abstract class BinaryFeatureMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> {

    public BinaryFeatureMapper() {

    }

    private static final int[] indices = new int[]{0, 0};

    @Override
    public int numberOfFeatures() {
        return Integer.bitCount(Integer.MAX_VALUE);
    }

    public abstract int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record);

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

    }

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
        int value = getIntegerValue(record);
        return (value & (1 << featureIndex)) != 0?1f:0f;
    }
}
