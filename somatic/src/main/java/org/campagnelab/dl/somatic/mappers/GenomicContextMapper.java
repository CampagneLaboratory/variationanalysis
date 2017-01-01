package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;
import java.util.function.Function;

/**
 * Maps the full genomic context using multiple onehotfeaturemapper
 * Created by rct66 on 10/25/16.
 */


public class GenomicContextMapper extends NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>, FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    private ConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    public GenomicContextMapper(Properties sbiProperties, int maxContextSize) {

        this(Math.min(maxContextSize, (int) Float.parseFloat(sbiProperties.getProperty("stats.genomicContextSize.min", "0.0"))));
        if (sbiProperties.getProperty("") == null) {
            throw new RuntimeException("Unable to obtain stats.genomicContextSize.min from properties.");
        }
    }

    public GenomicContextMapper(Properties sbiProperties) {

        this((int) Float.parseFloat(sbiProperties.getProperty("stats.genomicContextSize.min", "0.0")));
        if (sbiProperties.getProperty("stats.genomicContextSize.min") == null) {
            throw new RuntimeException("Unable to obtain stats.genomicContextSize.min from properties.");
        }
    }

    public GenomicContextMapper(int contextSize) {
        OneHotBaseFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>[] refContext = new OneHotBaseFeatureMapper[contextSize];
        for (int i = 0; i < contextSize; i++) {
            refContext[i] = new OneHotBaseFeatureMapper<>(i,
                    record -> trim(contextSize, record.getGenomicSequenceContext()));
        }
        delegate = new ConcatFeatureMapper<>(refContext);
    }

    /**
     * Trim a larger context to remain centered on the base of interest, but have only up to contextSize bases.
     *
     * @param contextSize                  target context size.
     * @param recordGenomicSequenceContext the genomic context from the sbi file, may be larger than contextSize
     * @return
     */
    String trim(int contextSize, String recordGenomicSequenceContext) {
       assert contextSize>=recordGenomicSequenceContext.length() :"The trim length must be smaller than the .sbi context length.";
        if (recordGenomicSequenceContext.length() == contextSize) {
            return recordGenomicSequenceContext;
        }
        int clipLength=(recordGenomicSequenceContext.length()-contextSize)/2;
        String result= recordGenomicSequenceContext.substring(clipLength,contextSize+clipLength);
        return result;
    }

    public GenomicContextMapper(int contextSize, Function<BaseInformationRecords.BaseInformationOrBuilder, String> function) {
        OneHotBaseFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>[] refContext = new OneHotBaseFeatureMapper[contextSize];
        for (int i = 0; i < contextSize; i++) {
            refContext[i] = new OneHotBaseFeatureMapper<>(i, function);
        }
        delegate = new ConcatFeatureMapper<>(refContext);
    }

    @Override
    public int numberOfFeatures() {
        return delegate.numberOfFeatures();

    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        delegate.prepareToNormalize(record, indexOfRecord);
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

    @Override
    public String getFeatureName(int featureIndex) {
        return "GenomicContextMapper" + featureIndex;
    }
}
