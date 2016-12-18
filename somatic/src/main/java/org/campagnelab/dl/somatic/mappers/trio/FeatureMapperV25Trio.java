package org.campagnelab.dl.somatic.mappers.trio;

import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.somatic.mappers.DensityMapper;
import org.campagnelab.dl.somatic.mappers.GenomicContextMapper;
import org.campagnelab.dl.somatic.mappers.NamingConcatFeatureMapper;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Properties;


/**
 * Same as V24, with no genomic position encoding.
 */
public class FeatureMapperV25Trio extends NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements ConfigurableFeatureMapper {
    private NamingConcatFeatureMapper delegate;

    private String recordTo(final int contextLength, BaseInformationRecords.BaseInformationOrBuilder record, int countIndex) {
        final List<BaseInformationRecords.CountInfo> counts = record.getSamples(record.getSamplesCount() - 1).getCountsList();
        ArrayList<BaseInformationRecords.CountInfo> sorted = new ArrayList<>();
        sorted.addAll(counts);
        Collections.sort(sorted, (o1, o2) ->
                (o2.getGenotypeCountForwardStrand() + o2.getGenotypeCountReverseStrand()) - (o1.getGenotypeCountForwardStrand() + o1.getGenotypeCountReverseStrand())
        );
        String to = countIndex < sorted.size() ? sorted.get(countIndex).getToSequence() : "N";
        StringBuffer context = new StringBuffer(to);

        while (context.length() < contextLength) {
            context.append("N");
        }
        context.setLength(contextLength);
        return context.toString();
    }

    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param sbiProperties properties from an sbi reader.
     */
    public void configure(Properties sbiProperties) {


        delegate = new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                new SimpleFeatureCalculatorTrio(true),
                new IndelFeaturesTrio(),
                new GenomicContextMapper(sbiProperties),
                // we reuse the genomic context mapper to map the to field for the first four genotypes (by decreasing count):
                new GenomicContextMapper(1, record -> recordTo(1, record, 0)),
                new GenomicContextMapper(1, record -> recordTo(1, record, 1)),
                new GenomicContextMapper(1, record -> recordTo(1, record, 2)),
                new GenomicContextMapper(1, record -> recordTo(1, record, 3)),
                new ReadIndexFeaturesTrio(),
                new FractionDifferences4Trio(0),
                new FractionDifferences4Trio(1),
                new MagnitudeFeatures2Trio(),
                new DensityMapper("numVariationsInRead", 20, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList)),
                new DensityMapper("readMappingQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                new DensityMapper("readMappingQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                new DensityMapper("baseQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                new DensityMapper("baseQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList)),
                new DensityMapper("insertSizes", 10, sbiProperties, (BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder) -> {
                    return TraversalHelper.forAllSampleCounts(baseInformationOrBuilder, BaseInformationRecords.CountInfo::getInsertSizesList);
                },
                        insertSize -> (float)Math.log10(insertSize))

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

}
