package org.campagnelab.dl.somatic.mappers;


import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.function.Function;

/**
 * Maps a int indexing into a record's genomic sequence context into a one hot base feature
 * Created by rct66 on 10/25/16.
 */
public class OneHotBaseMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
         {

    static private Logger LOG = LoggerFactory.getLogger(OneHotBaseMapper.class);

    private Function<BaseInformationRecords.BaseInformationOrBuilder, java.lang.String> recordToString;
    private int baseIndex;

    public OneHotBaseMapper(int baseIndex, Function<BaseInformationRecords.BaseInformationOrBuilder, java.lang.String> recordToString) {
        this.baseIndex = baseIndex;
        this.recordToString = recordToString;
    }

    private static final int[] indices = new int[]{0, 0};

    public int numberOfFeatures() {
        return 6;
    }

    public int getIntegerOfBase(BaseInformationRecords.BaseInformationOrBuilder record) {
        java.lang.String context = recordToString.apply(record);

        if (baseIndex < 0 || baseIndex >= context.length()) {
            LOG.warn("incompatible character index:" + baseIndex + " for context:" + context +" of length "+context.length());
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
        int value = getIntegerOfBase(record);
        return value == featureIndex ? 1F : 0F;
    }
}
