package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.datavec.api.records.reader.RecordReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Same as V18, but adds density features for numVariationsInRead. Starting to use Java8 lambdas to customize generic feature mappers.
 */
public class FeatureMapperV19 extends NamingConcatFeatureMapper implements ConfigurableFeatureMapper {
    NamingConcatFeatureMapper delegate;

    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param sbiProperties properties from an sbi reader.
     */
    public void configure(Properties sbiProperties) {

        delegate = new NamingConcatFeatureMapper(new SimpleFeatureCalculator(true),
                new IndelFeatures(),
                new ReadIndexFeaturesFix(),
                new FractionDifferences4(),
                new MagnitudeFeatures2(),
                new GenomicContextMapper(sbiProperties),
                new DensityMapper("numVariationsInRead", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)),
                new DensityMapper("readMappingQuality.forward", 25, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                new DensityMapper("readMappingQuality.reverse", 25, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                new DensityMapper("baseQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                new DensityMapper("baseQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList)));


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

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        delegate.mapFeatures(record, inputs, offset, indexOfRecord);
    }
}
