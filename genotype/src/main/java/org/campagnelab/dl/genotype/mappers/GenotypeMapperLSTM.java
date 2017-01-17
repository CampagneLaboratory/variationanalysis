package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Created by joshuacohen on 1/17/17.
 */
public class GenotypeMapperLSTM implements
        FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder>, ConfigurableFeatureMapper {
    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    @Override
    public void configure(Properties readerProperties) {
        delegate = new GenomicContextIndelMapper(readerProperties);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return delegate.getFeatureName(featureIndex);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();
    }

    @Override
    public MappedDimensions dimensions() {
        return delegate.dimensions();
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(record, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.isMasked(record, featureIndex);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }
}
