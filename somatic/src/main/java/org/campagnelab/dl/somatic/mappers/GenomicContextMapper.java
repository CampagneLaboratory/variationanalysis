package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.framework.iterators.ConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.somatic.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.somatic.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.somatic.utils.ProtoPredictor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Maps the full genomic context using multiple onehotfeaturemapper
 * Created by rct66 on 10/25/16.
 */


public class GenomicContextMapper extends AbstractFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> {
    private ConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    public GenomicContextMapper(Properties sbiProperties) {

        this((int) Float.parseFloat(sbiProperties.getProperty("stats.genomicContextSize.min", "0.0")));
        if (sbiProperties.getProperty("stats.genomicContextSize.min") == null) {
            throw new RuntimeException("Unable to obtain stats.genomicContextSize.min from properties.");
        }
    }

    public GenomicContextMapper(int contextSize) {
        OneHotBaseMapper[] refContext = new OneHotBaseMapper[contextSize];
        for (int i = 0; i < contextSize; i++) {
            refContext[i] = new OneHotBaseMapper(i, BaseInformationRecords.BaseInformationOrBuilder::getGenomicSequenceContext);
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
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        ReadIndexWithCounts myCounts = (ReadIndexWithCounts) count;
        myCounts.set(ProtoPredictor.expandFreq(sampleCounts.getReadIndicesForwardStrandList()),
                ProtoPredictor.expandFreq(sampleCounts.getReadIndicesReverseStrandList()));
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory() {
            @Override
            public GenotypeCount create() {
                return new ReadIndexWithCounts();
            }
        };
    }

}
