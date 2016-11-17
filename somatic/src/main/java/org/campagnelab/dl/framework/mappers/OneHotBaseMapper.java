package org.campagnelab.dl.framework.mappers;


import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Maps a int indexing into a record's genomic sequence context into a one hot base feature
 * Created by rct66 on 10/25/16.
 */
public class OneHotBaseMapper<RecordType> implements FeatureMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(OneHotBaseMapper.class);

    private Function<RecordType, String> recordToString;
    private int baseIndex;

    public OneHotBaseMapper(int baseIndex, Function<RecordType, String> recordToString) {
        this.baseIndex = baseIndex;
        this.recordToString = recordToString;
    }

    private static final int[] indices = new int[]{0, 0};

    public int numberOfFeatures() {
        return 6;
    }

    public int getIntegerOfBase(RecordType record) {
        String context = recordToString.apply(record);

        if (baseIndex < 0 || baseIndex >= context.length()) {
            LOG.warn("incompatible character index:" + baseIndex + " for context:" + context + " of length " + context.length());
            return 5;
        }
        Character base = context.charAt(baseIndex);
        int baseInt;
        switch (base) {
            case 'a':
            case 'A':
                baseInt = 0;
                break;
            case 't':
            case 'T':
                baseInt = 1;
                break;
            case 'c':
            case 'C':
                baseInt = 2;
                break;
            case 'g':
            case 'G':
                baseInt = 3;
                break;
            case 'n':
            case 'N':
                baseInt = 4;
                break;
            default:
                baseInt = 5;
                break;
        }
        return baseInt;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

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
        int value = getIntegerOfBase(record);
        return value == featureIndex ? 1F : 0F;
    }

    @Override
    public boolean hasMask() {
        return false;
    }

    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfFeatures());
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return false;
    }
}
