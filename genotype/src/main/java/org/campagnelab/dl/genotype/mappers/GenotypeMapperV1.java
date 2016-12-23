package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.MaxNormalizationMapper;
import org.campagnelab.dl.somatic.mappers.DensityMapper;
import org.campagnelab.dl.somatic.mappers.GenomicContextMapper;
import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * This m
 */
public class GenotypeMapperV1  extends GenotypeFeatureMapper{
    private NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;
    //default sampleIndex is zero, adjustable with setter
    private int sampleIndex = 0;


    public GenotypeMapperV1() {
        sortCounts = false;
    }

    /**
     * Configure the feature mapper to map a specific sampleIndex
     */
    public void configure(Properties sbiProperties) {

        FeatureNameMapper[] countMappers = new FeatureNameMapper[20];
        FeatureNameMapper[] readIndexMappers = new FeatureNameMapper[20];
        int genotypeIndex = 0;
        for (int i = 0; i < 10; i++) {
            countMappers[i] = (new SingleGenoTypeCountMapper(sampleIndex, i, true));
            readIndexMappers[i] = (new SingleReadIndexCountMapper(sampleIndex, i, true));
            genotypeIndex++;
        }
        genotypeIndex = 0;
        for (int i = 10; i < 20; i++) {

            countMappers[i] = (new SingleGenoTypeCountMapper(sampleIndex, genotypeIndex, false));
            readIndexMappers[i] = (new SingleReadIndexCountMapper(sampleIndex, genotypeIndex, false));
            genotypeIndex++;
        }
        delegate =
                new NamingConcatFeatureMapper<>(
                        new MaxNormalizationMapper(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(countMappers)),
                        new InverseNormalizationMapper(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(countMappers)),
                        new MaxNormalizationMapper(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readIndexMappers)),
                        new InverseNormalizationMapper(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readIndexMappers)),
                        new GenomicContextMapper(sbiProperties),
                        new DensityMapper("numVariationsInRead", 20, sbiProperties, baseInformationOrBuilder ->
                                TraversalHelper.forSampleCounts(sampleIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)),
                        new DensityMapper("readMappingQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                                TraversalHelper.forSampleCounts(sampleIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                        new DensityMapper("readMappingQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                                TraversalHelper.forSampleCounts(sampleIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                        new DensityMapper("baseQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                                TraversalHelper.forSampleCounts(sampleIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                        new DensityMapper("baseQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                                TraversalHelper.forSampleCounts(sampleIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList))
                )
        ;
        numFeatures = delegate.numberOfFeatures();

    }

    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    @Override
    public String getFeatureName(int i) {
        return delegate.getFeatureName(i);
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

}
