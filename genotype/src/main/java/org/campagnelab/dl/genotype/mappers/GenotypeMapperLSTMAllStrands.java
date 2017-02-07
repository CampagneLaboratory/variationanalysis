package org.campagnelab.dl.genotype.mappers;

import org.apache.commons.lang.math.NumberUtils;
import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.util.WarningCounter;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Properties;

/**
 * Created by joshuacohen on 2/6/17.
 */
public class GenotypeMapperLSTMAllStrands implements
        FeatureNameMapper<BaseInformationRecords.BaseInformation>, ConfigurableFeatureMapper {
    private static final Logger LOG = LoggerFactory.getLogger(GenotypeMapperLSTMAllStrands.class);
    private static final WarningCounter counter = new WarningCounter();
    private int sampleIndex;
    private int indelSequenceLength;
    private MappedDimensions dim;
    private int maskLen;
    private int[] indicesMapper = new int[]{0, 0, 0};
    private int[] indicesMasker = new int[]{0, 0};

    private Sample cachedSample = null;
    private static final int featuresPerSequence = 8;
    private static final int sequencesPerIndel = 4;
    private static final int featuresPerTimeStep = featuresPerSequence * sequencesPerIndel;
    private static final int defaultIndelSequenceLength = 30;

    @Override
    public void configure(Properties readerProperties) {
        String indelSequenceLengthProperty = readerProperties.getProperty("indelSequenceLength");
        if (indelSequenceLengthProperty == null) {
            indelSequenceLength = defaultIndelSequenceLength;
        } else {
            indelSequenceLength = Integer.parseInt(indelSequenceLengthProperty);
        }
        dim = new MappedDimensions(featuresPerTimeStep, indelSequenceLength);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return GenotypeMapperLSTMAllStrands.class.getSimpleName() + featureIndex;
    }

    @Override
    public int numberOfFeatures() {
        return dim.numElements();
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        BaseInformationRecords.SampleInfo sampleInfo = record.getSamples(sampleIndex);
        cachedSample = new Sample(sampleInfo.getCounts(0).getFromSequence(),
                sampleInfo.getCounts(0).getToSequence(),
                sampleInfo.getCounts(1).getToSequence(),
                sampleInfo.getCounts(2).getToSequence(),
                indelSequenceLength);
        maskLen = cachedSample.maxLen();
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformation record, INDArray inputs, int indexOfRecord) {
        indicesMapper[0] = indexOfRecord;
        for (int i = 0; i < indelSequenceLength; i++) {
            indicesMapper[2] = i;
            for (int j = 0; j < featuresPerTimeStep; j++) {
                indicesMapper[1] = j;
                int featureIndex = i * featuresPerTimeStep + j;
                inputs.putScalar(indicesMapper, produceFeature(record, featureIndex));
            }
        }
    }

    @Override
    public boolean hasMask() {
        return true;
    }

    @Override
    public void maskFeatures(BaseInformationRecords.BaseInformation record, INDArray mask, int indexOfRecord) {
        indicesMasker[0] = indexOfRecord;
        for (int i = 0; i < indelSequenceLength; i++) {
            indicesMasker[1] = i;
            for (int j = 0; j < featuresPerTimeStep; j++) {
                int featureIndex = i * featuresPerTimeStep + j;
                mask.putScalar(indicesMasker, isMasked(record, featureIndex) ? 1F : 0F);
            }
        }
    }

    @Override
    public boolean isMasked(BaseInformationRecords.BaseInformation record, int featureIndex) {
        return featureIndex / featuresPerTimeStep < maskLen;
    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformation record, int featureIndex) {
        int featureInTimeStepIndex = featureIndex % featuresPerTimeStep;
        int sequenceIndex = featureInTimeStepIndex / featuresPerSequence;
        int featureInSequenceIndex = featureInTimeStepIndex % featuresPerSequence;
        String sequence = cachedSample.getSequence(sequenceIndex);
        int timeStepIndex = featureIndex / featuresPerTimeStep;
        if (timeStepIndex < sequence.length()) {
            int sequenceBase = getIntegerOfBase(sequence, timeStepIndex);
            return sequenceBase == featureInSequenceIndex ? 1F : 0F;
        } else if (timeStepIndex < maskLen) {
            return featureInSequenceIndex == featuresPerSequence - 2 ? 1F : 0F;
        } else {
            return 0f;
        }
    }

    private static int getIntegerOfBase(String field, int baseIndex) {
        if (baseIndex < 0 || baseIndex >= field.length()) {
            counter.warn(LOG, String.format("incompatible character index: %c for field: %s of length %d",
                    baseIndex, field, field.length()));
            return 5;
        }
        Character base = field.charAt(baseIndex);
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
                baseInt = 7;
                break;
        }
        return baseInt;
    }


    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    private class Sample {
        private final int maxLen;
        private final String from;
        private final String to1;
        private final String to2;
        private final String to3;

        private Sample(String from, String to1, String to2, String to3, int maxLen) {
            this.from = from;
            this.to1 = to1;
            this.to2 = to2;
            this.to3 = to3;
            this.maxLen = maxLen;
        }

        private String getSequence(int sequenceIndex) {
            switch (sequenceIndex) {
                case 0:
                    return from;
                case 1:
                    return to1;
                case 2:
                    return to2;
                case 3:
                    return to3;
                default:
                    throw new IllegalArgumentException("Invalid sequence idx " + sequenceIndex);
            }
        }

        private int maxLen() {
            int seqMaxLen = NumberUtils.max(new int[]{from.length(), to1.length(), to2.length(), to3.length()});
            return Math.min(maxLen, seqMaxLen);
        }
    }
}
