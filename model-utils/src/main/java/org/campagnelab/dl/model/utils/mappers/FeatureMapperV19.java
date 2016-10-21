package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
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
     * @param reader
     */
    public void configure(SequenceBaseInformationReader reader) {
        Properties sbiProperties = reader.getProperties();
        if (!propertiesPresent(sbiProperties, "stats.numVariationsInRead")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for numVariationsInRead (stats.numVariationsInRead.min and stats.numVariationsInRead.max)");
        }
        float minNumVariationsInRead = getMin(sbiProperties, "stats.numVariationsInRead");
        float maxNumVariationsInRead = getMax(sbiProperties, "stats.numVariationsInRead");

        if (!propertiesPresent(sbiProperties, "stats.readMappingQuality.forward")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for stats.baseQuality.forward");
        }
        float minMappingQualityForward = getMin(sbiProperties, "stats.readMappingQuality.forward");
        float maxMappingQualityForward = getMax(sbiProperties, "stats.readMappingQuality.forward");

        if (!propertiesPresent(sbiProperties, "stats.readMappingQuality.reverse")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for stats.baseQuality.forward");
        }
        float minMappingQualityReverse = getMin(sbiProperties, "stats.readMappingQuality.reverse");
        float maxMappingQualityReverse = getMax(sbiProperties, "stats.readMappingQuality.reverse");

        if (!propertiesPresent(sbiProperties, "stats.baseQuality.forward")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for stats.baseQuality.forward");
        }
        float minBaseQualityForward = getMin(sbiProperties, "stats.baseQuality.forward");
        float maxBaseQualityForward = getMax(sbiProperties, "stats.baseQuality.forward");

        if (!propertiesPresent(sbiProperties, "stats.baseQuality.reverse")) {
            throw new UnsupportedOperationException("The sbip file does not contain the statistics for stats.baseQuality.forward");
        }
        float minBaseQualityReverse = getMin(sbiProperties, "stats.baseQuality.reverse");
        float maxBaseQualityReverse = getMax(sbiProperties, "stats.baseQuality.reverse");

        delegate = new NamingConcatFeatureMapper(new SimpleFeatureCalculator(true),
                new IndelFeatures(),
                new ReadIndexFeaturesFix(),
                new FractionDifferences4(),
                new MagnitudeFeatures2(),
                new DensityMapper("numVariationsInRead", 10, minNumVariationsInRead, maxNumVariationsInRead, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)),
                new DensityMapper("mappingQualityForward", 25, minMappingQualityForward, maxMappingQualityForward, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                new DensityMapper("mappingQualityReverse", 25, minMappingQualityReverse, maxMappingQualityReverse, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                new DensityMapper("baseQualityForward", 10, minBaseQualityForward, maxBaseQualityForward, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                new DensityMapper("baseQualityReverse", 10, minBaseQualityReverse, maxBaseQualityReverse, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList)));

    }

    private float getMin(Properties sbiProperties, String propertyName) {
        return Float.parseFloat(sbiProperties.getProperty(propertyName + ".min"));
    }

    private float getMax(Properties sbiProperties, String propertyName) {
        return Float.parseFloat(sbiProperties.getProperty(propertyName + ".max"));
    }

    private boolean propertiesPresent(Properties sbiProperties, String s) {
        return sbiProperties.containsKey(s + ".min") || !sbiProperties.containsKey(s + ".max");
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
