package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Maps an integer to 32 binary features, using one hot encoding.
 * Created by fac2003 on 7/12/16.
 */
public abstract class BinaryFeatureMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    private int maxValue;
    private int value;

    public BinaryFeatureMapper(int maxValue) {
        this.maxValue = maxValue;
    }

    public BinaryFeatureMapper() {
        this(Integer.MAX_VALUE);
    }

    private static final int[] indices = new int[]{0, 0};

    @Override
    public int numberOfFeatures() {
        return Integer.bitCount(maxValue);
    }

    public abstract int getIntegerValue(BaseInformationRecords.BaseInformationOrBuilder record);

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        value = getIntegerValue(record);
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

        return (value >> featureIndex & 1);
    }
}
