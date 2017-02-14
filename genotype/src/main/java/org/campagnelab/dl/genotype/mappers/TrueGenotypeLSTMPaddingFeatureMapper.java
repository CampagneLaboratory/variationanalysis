package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;
import java.util.function.Function;

/**
 * Created by joshuacohen on 2/8/17.
 */
public class TrueGenotypeLSTMPaddingFeatureMapper implements
        FeatureNameMapper<BaseInformationRecords.BaseInformation>, ConfigurableFeatureMapper {

    private RNNFeatureMapper<String> delegate;
    private String cachedRecordGenotype;
    private int trueGenotypeLength;
    private static final int defaultIndelSequenceLength = 30;

    @Override
    public void configure(Properties readerProperties) {
        String trueGenotypeLengthProperty = readerProperties.getProperty("trueGenotypeLength");
        if (trueGenotypeLengthProperty == null) {
            trueGenotypeLength = defaultIndelSequenceLength;
        } else {
            trueGenotypeLength = Integer.parseInt(trueGenotypeLengthProperty);
        }
        OneHotBaseFeatureMapper<String>[] delegateMapperArray = new OneHotBaseFeatureMapper[trueGenotypeLength];
        for (int i = 0; i < trueGenotypeLength; i++) {
            delegateMapperArray[i] = new OneHotBaseFeatureMapper<>(i, Function.identity(),
                   TrueGenotypeLSTMPaddingFeatureMapper::baseToPaddingFeature, 1);
        }
        delegate = new RNNFeatureMapper<>(String::length, delegateMapperArray);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return TrueGenotypeLSTMPaddingFeatureMapper.class.getSimpleName() + featureIndex;
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
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        cachedRecordGenotype = record.getTrueGenotype();
        delegate.prepareToNormalize(cachedRecordGenotype, indexOfRecord);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformation record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(cachedRecordGenotype, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(BaseInformationRecords.BaseInformation record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(cachedRecordGenotype, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(BaseInformationRecords.BaseInformation record, int featureIndex) {
        return delegate.isMasked(cachedRecordGenotype, featureIndex);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformation record, int featureIndex) {
        return delegate.produceFeature(cachedRecordGenotype, featureIndex);
    }

    private static int baseToPaddingFeature(String recordString, int baseIndex) {
        return 0;
    }
}
