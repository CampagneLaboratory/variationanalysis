package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.util.WarningCounter;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Properties;
import java.util.function.Function;

/**
 * Created by joshuacohen on 1/17/17.
 */
public class GenotypeMapperLSTM implements
        FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder>, ConfigurableFeatureMapper {
    private int sampleIndex;
    private RNNFeatureMapper<String> delegate;
    private String cachedRecordIndelString;
    private int indelSequenceLength;

    public enum Input {
        FROM,
        G1,
        G2,
        G3,
    }
    private Input inputType = Input.FROM;
    private static final int featuresPerOHBM = 7;

    @Override
    public void configure(Properties readerProperties) {
        indelSequenceLength = Integer.parseInt(readerProperties.getProperty("indelSequenceLength"));
        OneHotBaseFeatureMapper<String>[] delegateMapperArray = new OneHotBaseFeatureMapper[indelSequenceLength];
        for (int i = 0; i < indelSequenceLength; i++) {
            delegateMapperArray[i] = new OneHotBaseFeatureMapper<>(i, Function.identity(),
                    GenotypeMapperLSTM::getIntegerOfBase, featuresPerOHBM);
        }
        delegate = new RNNFeatureMapper<>(String::length, delegateMapperArray);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return GenotypeMapperLSTM.class.getSimpleName() + featureIndex;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        BaseInformationRecords.SampleInfo sampleInfo = record.getSamples(sampleIndex);
        switch (inputType) {
            case FROM:
                cachedRecordIndelString = sampleInfo.getCounts(0).getFromSequence();
                break;
            case G1:
                cachedRecordIndelString = sampleInfo.getCounts(0).getToSequence();
                break;
            case G2:
                cachedRecordIndelString = sampleInfo.getCounts(1).getToSequence();
                break;
            case G3:
                cachedRecordIndelString = sampleInfo.getCounts(2).getToSequence();
                break;
            default:
                throw new RuntimeException("Invalid input type");
        }
        if (cachedRecordIndelString.length() > indelSequenceLength) {
            cachedRecordIndelString = cachedRecordIndelString.substring(0, indelSequenceLength);
        }
        delegate.prepareToNormalize(cachedRecordIndelString, indexOfRecord);
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
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        delegate.mapFeatures(cachedRecordIndelString, inputs, indexOfRecord);
    }

    @Override
    public boolean hasMask() {
        return delegate.hasMask();
    }

    @Override
    public void maskFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray mask, int indexOfRecord) {
        delegate.maskFeatures(cachedRecordIndelString, mask, indexOfRecord);
    }

    @Override
    public boolean isMasked(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.isMasked(cachedRecordIndelString, featureIndex);
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(cachedRecordIndelString, featureIndex);
    }

    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }


    public void setInputType(Input inputType) {
        this.inputType = inputType;
    }

    private static WarningCounter counter = new WarningCounter();
    private static final Logger LOG = LoggerFactory.getLogger(GenotypeMapperLSTM.class);

    private static int getIntegerOfBase(String context, int baseIndex) {
        if (baseIndex < 0 || baseIndex >= context.length()) {
            counter.warn(LOG, "incompatible character index: %d for context: {} of length {}",
                    baseIndex, context, context.length());
            return 6;
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
            case '-':
                baseInt = 5;
                break;
            default:
                baseInt = 6;
                break;
        }
        return baseInt;
    }
}
