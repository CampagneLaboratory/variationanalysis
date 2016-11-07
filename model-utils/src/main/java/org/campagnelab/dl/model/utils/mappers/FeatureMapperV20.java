package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Same as V19, but more point for some density mappers (numVariationsInRead).
 */
public class FeatureMapperV20 extends NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements ConfigurableFeatureMapper {
    private NamingConcatFeatureMapper delegate;

    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param sbiProperties properties from an sbi reader.
     */
    public void configure(Properties sbiProperties) {

        delegate = new NamingConcatFeatureMapper(
                new SimpleFeatureCalculator(true),
                new IndelFeatures(),
                new GenomicContextMapper(sbiProperties),
                new ReadIndexFeaturesFix(),
                new FractionDifferences4(),
                new MagnitudeFeatures2(),
                new DensityMapper("numVariationsInRead", 20, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)),
                new DensityMapper("readMappingQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                new DensityMapper("readMappingQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                new DensityMapper("baseQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                new DensityMapper("baseQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList))
        );


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
