package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;
import java.util.function.Function;

/**
 * Created by joshuacohen on 2/8/17.
 */
public class TrueGenotypeLSTMDecodingFeatureMapper implements
        FeatureNameMapper<BaseInformationRecords.BaseInformation>, ConfigurableFeatureMapper {

    private RNNFeatureMapper<String> delegate;
    private String cachedRecordGenotype;
    private int trueGenotypeLength;
    private static final int defaultGenotypeSequenceLength = 30;
    private boolean isPredicting;

    @Override
    public void configure(Properties readerProperties) {
        isPredicting = Boolean.parseBoolean(readerProperties.getProperty("isPredicting"));
        String trueGenotypeLengthProperty = readerProperties.getProperty("trueGenotypeLength");
        if (trueGenotypeLengthProperty == null) {
            trueGenotypeLength = defaultGenotypeSequenceLength;
        } else {
            trueGenotypeLength = Integer.parseInt(trueGenotypeLengthProperty);
        }
        OneHotBaseFeatureMapper<String>[] delegateMapperArray = new OneHotBaseFeatureMapper[trueGenotypeLength + 3];
        for (int i = 0; i < trueGenotypeLength + 3; i++) {
            delegateMapperArray[i] = new OneHotBaseFeatureMapper<>(i, Function.identity(),
                   TrueGenotypeLSTMDecodingFeatureMapper::baseToPaddingFeature,
                    TrueGenotypeLSTMLabelMapper.featuresOrLabelsPerTimeStep);
        }
        delegate = new RNNFeatureMapper<>(String::length, delegateMapperArray);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return TrueGenotypeLSTMDecodingFeatureMapper.class.getSimpleName() + featureIndex;
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
        String trueGenotype = record.getTrueGenotype();
        StringBuilder cachedRecordGenotypeBuilder = new StringBuilder();
        cachedRecordGenotypeBuilder.append("$$");
        if (!isPredicting) {
            if (trueGenotype.length() >= trueGenotypeLength) {
                cachedRecordGenotypeBuilder.append(trueGenotype.substring(0, trueGenotypeLength));
            } else {
                cachedRecordGenotypeBuilder.append(trueGenotype);
            }
            cachedRecordGenotypeBuilder.append('*');
        }
        cachedRecordGenotype = cachedRecordGenotypeBuilder.toString();
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
        return TrueGenotypeLSTMLabelMapper.baseToLabel(recordString.charAt(baseIndex));
    }
}
