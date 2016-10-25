package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;

/**
 * Maps the full genomic context using multiple onehotfeaturemapper
 * Created by rct66 on 10/25/16.
 */



public class GenomicContextMapper implements FeatureMapper {
    private OneHotBaseMapper[] refContext;
    private ConcatFeatureMapper delegate;

    public GenomicContextMapper(int contextSize) {
        refContext = new OneHotBaseMapper[contextSize];
        for (int i = 0; i < contextSize; i++){
            refContext[i] = new OneHotBaseMapper(i);
        }
        delegate = new ConcatFeatureMapper(refContext);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }
}
