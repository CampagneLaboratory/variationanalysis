package org.campagnelab.dl.framework.mappers;


import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Maps an object to a fixed set of bit elements. Transform the object to a hashcode and reduce with modulo
 * to a fixed number of values. Then use one hot encoding to yield a vector.
 */
public class OneHotHashModuloMapper<RecordType> extends NoMaskFeatureMapper<RecordType> {

    static private Logger LOG = LoggerFactory.getLogger(OneHotHashModuloMapper.class);

    private Function<RecordType, Object> recordToString;
    // the number of elements in the output vector.
    private int vectorNumElements;
    // populated by prepareToNormalize
    private int reducedValue = -1;

    public OneHotHashModuloMapper(int vectorNumElements, Function<RecordType, Object> recordToString) {
        this.vectorNumElements = vectorNumElements;
        this.recordToString = recordToString;
    }

    private static final int[] indices = new int[]{0, 0};

    public int numberOfFeatures() {
        return vectorNumElements;
    }

    public int getIntegerOfBase(RecordType record) {
        Object context = recordToString.apply(record);
        return context.hashCode() % vectorNumElements;

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
        assert reducedValue >= 0: "prepareToNormalize must be called before produceFeature.";
        return reducedValue == featureIndex ? 1F : 0F;
    }
}
