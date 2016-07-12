package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A feature mapper that reduces an integer to a cateogy and encodes with one hot encoding.
 * Created by fac2003 on 7/12/16.
 */
public abstract class CategoryFeatureMapper implements FeatureMapper {
    private int minCategoryIndex;
    private int maxCategoryIndex;

    public CategoryFeatureMapper(int minCategoryIndex, int maxCategoryIndex) {
        assert minCategoryIndex>0;
        assert maxCategoryIndex>=minCategoryIndex;
        this.minCategoryIndex = minCategoryIndex;
        this.maxCategoryIndex = maxCategoryIndex;
    }

    private static final int[] indices = new int[]{0, 0};

    @Override
    public int numberOfFeatures() {
        return maxCategoryIndex - minCategoryIndex + 1;
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
        assert value<=maxCategoryIndex : String.format("value %d cannot be more than maxCategoryIndex(%d)", value,maxCategoryIndex);
        value=Math.min(value,maxCategoryIndex);
        return value-minCategoryIndex;
    }
}
