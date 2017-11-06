package org.campagnelab.dl.framework.mappers;


import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Maps an integer to a fixed set of values. Uses one hot encoding to yield a vector.
 */
public class OneHotIntegerMapper<RecordType> extends NoMaskFeatureMapper<RecordType> {

    static private Logger LOG = LoggerFactory.getLogger(OneHotIntegerMapper.class);

    private Function<RecordType, Integer> recordToInteger;
    // the number of elements in the output vector.
    private int vectorNumElements;
    // populated by prepareToNormalize
    private int reducedValue = -1;

    public OneHotIntegerMapper(int vectorNumElements, Function<RecordType, Integer> recordToString) {
        this.vectorNumElements = vectorNumElements;
        this.recordToInteger = recordToString;
    }

    private static final int[] indices = new int[]{0, 0};

    public int numberOfFeatures() {
        return vectorNumElements;
    }

    public int getIntegerOfBase(RecordType record) {
        int value = recordToInteger.apply(record);
        assert value < vectorNumElements : String.format("value %d cannot exceed vectorNumElements %d", value, vectorNumElements);
        return value;

    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {
        reducedValue = getIntegerOfBase(record);
    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {

        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        assert reducedValue >= 0 : "prepareToNormalize must be called before produceFeature.";
        return reducedValue == featureIndex ? 1F : 0F;
    }
}
