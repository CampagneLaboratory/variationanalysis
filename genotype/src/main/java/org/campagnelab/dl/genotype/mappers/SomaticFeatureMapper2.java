package org.campagnelab.dl.genotype.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.somatic.mappers.*;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;
import java.util.Set;


/**
 * A somatic feature mapper that reuses GenotypeMapperV37 and adds somatic specific features.
 */
public class SomaticFeatureMapper2 extends NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>
        implements ConfigurableFeatureMapper {
    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;

    private String recordTo(final int contextLength, BaseInformationRecords.BaseInformationOrBuilder record, int countIndex) {
        return MappingFunctions.recordTo(contextLength, record, countIndex);
    }

    int MAX_GENOTYPES;
    boolean withCombinedLayer;
    boolean withDistinctAlleleCounts;
    int germlineIndex;
    int somaticIndex;


    public SomaticFeatureMapper2() {
        this(0, 1);
        System.out.println("Somatic Feature Mapper instantiated with defaults:\n" +
                "germline = sample 0 in sbi/protobuf, somatic = sample 1 in sbi/protobuf");
    }

    public SomaticFeatureMapper2(int germlineIndex, int somaticIndex) {
        super();
        this.germlineIndex = germlineIndex;
        this.somaticIndex = somaticIndex;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        MAX_GENOTYPES = 3;
    }


    /**
     * Configure the feature mapper for a specific set of sbi files. This method accesses the properties of the reader.
     *
     * @param sbiProperties properties from an sbi reader.
     */
    public void configure(Properties sbiProperties) {
        String ploidyString = sbiProperties.getProperty("genotypes.ploidy");
        if (ploidyString == null) {
            MAX_GENOTYPES = 4;
        } else {
            MAX_GENOTYPES = Integer.parseInt(ploidyString) + 2;
        }

        String genomicContextLengthString = sbiProperties.getProperty("stats.genomicContextSize.min");
        assert genomicContextLengthString != null : "property must exist: stats.genomicContextSize.min";
        int genomicContextLength = (int) Float.parseFloat(genomicContextLengthString);

        String indelSequenceLengthString = sbiProperties.getProperty("indelSequenceLength");
        assert genomicContextLengthString != null : "property must exist: indelSequenceLength";
        int indelSequenceLength = Integer.parseInt(indelSequenceLengthString);
        final int indelMappedLength = indelSequenceLength;

        Set<Integer> sampleIndices = new ObjectArraySet<>();
        sampleIndices.add(germlineIndex);
        sampleIndices.add(somaticIndex);

        FeatureNameMapper[] matchesRefMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] firstBaseMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] originalGobyCountIndexMappers = new FeatureNameMapper[MAX_GENOTYPES];

        int genotypeIndex = 0;

        for (int i = 0; i < MAX_GENOTYPES; i++) {
            final int constantGenotypeIndex = genotypeIndex;

            matchesRefMappers[i] = (new MatchesReferenceMapper(somaticIndex, i));
            firstBaseMappers[i] = new GenomicContextMapper(indelMappedLength,
                    record -> {
                        final String toSequence = record.getSamples(0).getCounts(constantGenotypeIndex).getToSequence();
                        return toSequence.substring(0, Math.min(indelMappedLength, toSequence.length()));
                    }, true /* no warning if index outside of context, needed since indels have variable lengths */);
            originalGobyCountIndexMappers[i] = new OriginalGobyCountIndexMapper(somaticIndex, constantGenotypeIndex);
        }

        final OneSampleMapperUnsortedV2 a = new OneSampleMapperUnsortedV2(germlineIndex);
        final OneSampleMapperUnsortedV2 b = new OneSampleMapperUnsortedV2(somaticIndex);
        a.configure(sbiProperties);
        b.configure(sbiProperties);
        delegate = new CountReorderingMapper(somaticIndex, new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                a,
                b,
                /** Map the from sequence for genotypeIndex=0: */
                new GenomicContextMapper(indelMappedLength,
                        record -> {
                            final String fromSequence = record.getSamples(0).getCounts(0).getFromSequence();
                            return fromSequence.substring(0, Math.min(indelMappedLength, fromSequence.length()));
                        }, true /* no warning if index outside of context, needed since indels have variable lengths */), //same
                new GenomicContextMapper(sbiProperties),
                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(matchesRefMappers),
                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(originalGobyCountIndexMappers),
                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(firstBaseMappers),
         /* NumVariationsInReads for counts not in the best 3: */
                new DensityMapper("numVariationsInRead",
                        10, sbiProperties,
                        record -> TraversalHelper.forAllSampleCounts(record,
                                BaseInformationRecords.CountInfoOrBuilder::getNumVariationsInReadsList)),

                new DensityMapper("readMappingQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forNSampleCounts(sampleIndices, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList)),
                new DensityMapper("readMappingQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forNSampleCounts(sampleIndices, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList)),
                new DensityMapper("baseQuality.forward", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forNSampleCounts(sampleIndices, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList)),
                new DensityMapper("baseQuality.reverse", 10, sbiProperties, baseInformationOrBuilder ->
                        TraversalHelper.forNSampleCounts(sampleIndices, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList)),
                new DensityMapper("insertSizes", 10, sbiProperties, (BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder) -> {
                    return TraversalHelper.forNSampleCounts(sampleIndices, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getInsertSizesList);
                },
                        insertSize -> (float) Math.log10(insertSize)),
                new FractionDifferences4(),
                new MagnitudeFeatures2()
        ));

        numFeatures = delegate.numberOfFeatures();

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
        assert record.getSamplesCount() >= 2 : "A minimum of two samples is needed in the .sbi file for the somatic mapper.";
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
