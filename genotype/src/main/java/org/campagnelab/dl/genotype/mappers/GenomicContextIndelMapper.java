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
 * Maps the full genomic context using multiple onehotfeaturemapper
 * Created by rct66 on 10/25/16.
 */


public class GenomicContextIndelMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>, FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    private RNNFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    public GenomicContextIndelMapper(Properties sbiProperties, int maxContextSize) {

        this(Math.min(maxContextSize, (int) Float.parseFloat(sbiProperties.getProperty("stats.genomicContextSize.min", "0.0"))));
        if (sbiProperties.getProperty("") == null) {
            throw new RuntimeException("Unable to obtain stats.genomicContextSize.min from properties.");
        }
    }

    public GenomicContextIndelMapper(Properties sbiProperties) {
        this((int) Float.parseFloat(sbiProperties.getProperty("stats.genomicContextSize.min", "0.0")));
        if (sbiProperties.getProperty("stats.genomicContextSize.min") == null) {
            throw new RuntimeException("Unable to obtain stats.genomicContextSize.min from properties.");
        }
    }


    public GenomicContextIndelMapper(int contextSize) {
        OneHotBaseFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>[] refContext = new OneHotBaseFeatureMapper[contextSize];
        for (int i = 0; i < contextSize; i++) {
            refContext[i] = new OneHotBaseFeatureMapper<>(i,
                    record -> trim(contextSize, record.getGenomicSequenceContext()),
                    GenomicContextIndelMapper::getIntegerOfBase, 7);
        }
        delegate = new RNNFeatureMapper<>(r -> contextSize, refContext);
    }

    /**
     * Trim a larger context to remain centered on the base of interest, but have only up to trimLength bases.
     *
     * @param trimLength                  target context size.
     * @param recordGenomicSequenceContext the genomic context from the sbi file, may be larger than trimLength
     * @return
     */
    String trim(int trimLength, String recordGenomicSequenceContext) {
        assert trimLength<=recordGenomicSequenceContext.length() :
                String.format("The trim length (%d) must be smaller than the .sbi context length (%d).",
                        trimLength,  recordGenomicSequenceContext.length()) ;
        if (recordGenomicSequenceContext.length() == trimLength) {
            return recordGenomicSequenceContext;
        }
        int clipLength=(recordGenomicSequenceContext.length()-trimLength)/2;
        return recordGenomicSequenceContext.substring(clipLength,trimLength+clipLength);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();

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
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return delegate.produceFeature(record, featureIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        return "GenomicContextIndelMapper" + featureIndex;
    }

    private static WarningCounter counter = new WarningCounter();
    private static final Logger LOG = LoggerFactory.getLogger(GenomicContextIndelMapper.class);

    private static int getIntegerOfBase(String context, int baseIndex) {
        if (baseIndex < 0 || baseIndex >= context.length()) {
            counter.warn(LOG, String.format("incompatible character index: {} for context: {} of length {}",
                    baseIndex, context, context.length()));
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

