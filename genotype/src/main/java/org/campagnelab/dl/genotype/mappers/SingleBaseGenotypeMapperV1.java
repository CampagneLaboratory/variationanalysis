package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.somatic.mappers.*;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.CountInfoOrBuilder;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * Based on V37, but restricted to a single base of genotype and indel context.
 */
public class SingleBaseGenotypeMapperV1 extends GenotypeMapperV11 {


    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;
    //default sampleIndex is zero, adjustable with setter
    private int sampleIndex = 0;
    private OneHotBaseFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> offsetMapper;

    public SingleBaseGenotypeMapperV1() {
        this(0);
    }

    public SingleBaseGenotypeMapperV1(int sampleIndex) {
        super();
        this.sampleIndex = sampleIndex;
        sortCounts = true;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        MAX_GENOTYPES = 3;
        withSoftmaxGenotype = true;
    }

    int offset;

    public void setOffset(int offset) {
        this.offset = offset;
    }

    /**
     * Configure the feature mapper to map a specific sampleIndex
     */
    public void configure(Properties sbiProperties) {
        String ploidyString = sbiProperties.getProperty("genotypes.ploidy");
        if (ploidyString == null) {
            MAX_GENOTYPES = 3;
        } else {
            MAX_GENOTYPES = Integer.parseInt(ploidyString) + 1;
        }

        String genomicContextLengthString = sbiProperties.getProperty("stats.genomicContextSize.min");
        assert genomicContextLengthString != null : "property must exist: stats.genomicContextSize.min";
        int genomicContextLength = (int) Float.parseFloat(genomicContextLengthString);
        assert genomicContextLength == 1 : "Single base genotype mapper requires genomicContextLength property=1";

        String indelSequenceLengthString = sbiProperties.getProperty("indelSequenceLength");
        assert indelSequenceLengthString != null : "property must exist: indelSequenceLength";
        int indelSequenceLength = Integer.parseInt(indelSequenceLengthString);
        assert indelSequenceLength == 1 : "Single base genotype mapper requires indelSequenceLength property=1";

        final int indelMappedLength = indelSequenceLength;
        FeatureNameMapper[] matchesRefMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] countMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] readIndexMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] readMappingQualityMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] baseQualityMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] firstBaseMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] numVariationsInReadMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] targetAlignedLengthMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] queryAlignedLengthMappers = new FeatureNameMapper[MAX_GENOTYPES];

        //combined distances for now
        FeatureNameMapper[] distancesToReadVariations = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] distancesFromStartOfRead = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] distancesFromEndOfRead = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] bamFlagMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] originalGobyCountIndexMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] queryPositions = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] countOffsets = new FeatureNameMapper[MAX_GENOTYPES];


        int genotypeIndex = 0;

        for (int i = 0; i < MAX_GENOTYPES; i++) {
            final int constantGenotypeIndex = genotypeIndex;
            countMappers[i] = (new SingleGenoTypeCountMapper(sampleIndex, i, true));
            readIndexMappers[i] = (new SingleReadIndexCountMapper(sampleIndex, i, true));
            matchesRefMappers[i] = (new MatchesReferenceMapper(sampleIndex, i));
            FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> oneHotOffsetMapper =
                    new OneHotIntegerMapper<>(50,
                            record -> record.getSamples(sampleIndex).getCounts(constantGenotypeIndex).getOffset()
                    );
            countOffsets[i] =
                    new NamedWrapper<BaseInformationRecords.BaseInformationOrBuilder>(oneHotOffsetMapper) {

                        @Override
                        public String getFeatureName(int featureIndex) {
                            return "offset";
                        }
                    };


            firstBaseMappers[i] = new GenomicContextMapper(indelMappedLength,
                    record -> {
                        final String toSequence = record.getSamples(sampleIndex).getCounts(constantGenotypeIndex).getToSequence();
                        return toSequence.substring(0, Math.min(indelMappedLength, toSequence.length()));
                    }, true /* no warning if index outside of context, needed since indels have variable lengths */);

            queryPositions[i] = new DensityMapper("queryPosition",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQueryPositionsList) /*,
                    queryPosition -> (float)(Math.log(queryPosition+1)/Math.log(2))*/);

            numVariationsInReadMappers[i] = new DensityMapper("numVariationsInRead",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList));
            //bin width 1 density mapper that ignores variations outside of caps
            distancesToReadVariations[i] = new DensityMapperCapped("distancesToReadVariations.forward", "distancesToReadVariations.reverse",
                    -50, 50, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotypeBothStrands(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder,
                                    BaseInformationRecords.CountInfo::getDistancesToReadVariationsForwardStrandList,
                                    BaseInformationRecords.CountInfo::getDistancesToReadVariationsReverseStrandList));

            //bin width 1 density mapper that focuses on half the size of the genomic context (only length where we can observe homopolymers):
            distancesFromStartOfRead[i] = new DensityMapperCapped("distanceToStartOfRead",
                    0, genomicContextLength / 2, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder,
                                    BaseInformationRecords.CountInfo::getDistanceToStartOfReadList));

            //bin width 1 density mapper that focuses on half the size of the genomic context (only length where we can observe homopolymers):
            distancesFromEndOfRead[i] = new DensityMapperCapped("distanceToEndOfRead",
                    0, genomicContextLength / 2, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder,
                                    BaseInformationRecords.CountInfo::getDistanceToEndOfReadList));

            readMappingQualityMappers[i] = new DensityMapper("readMappingQuality.forward",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex,
                                    baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityForwardStrandList));
            baseQualityMappers[i] = new DensityMapper("baseQuality.forward",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresForwardStrandList));

            targetAlignedLengthMappers[i] = new DensityMapper("targetAlignedLength",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getTargetAlignedLengthsList));
            queryAlignedLengthMappers[i] = new DensityMapper("queryAlignedLength",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQueryAlignedLengthsList));
            bamFlagMappers[i] = new BamFlagMapper(sampleIndex, genotypeIndex);
            originalGobyCountIndexMappers[i] = new OriginalGobyCountIndexMapper(sampleIndex, constantGenotypeIndex);
            genotypeIndex++;
        }
        genotypeIndex = 0;
        for (int i = MAX_GENOTYPES; i < MAX_GENOTYPES * 2; i++) {
            final int constantGenotypeIndex = genotypeIndex;

            countMappers[i] = (new SingleGenoTypeCountMapper(sampleIndex, genotypeIndex, false));
            readIndexMappers[i] = (new SingleReadIndexCountMapper(sampleIndex, genotypeIndex, false));

            readMappingQualityMappers[i] = new DensityMapper("readMappingQuality.reverse",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex,
                                    baseInformationOrBuilder, BaseInformationRecords.CountInfo::getReadMappingQualityReverseStrandList));
            baseQualityMappers[i] = new DensityMapper("baseQuality.reverse",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQualityScoresReverseStrandList));
            genotypeIndex++;
        }


        delegate =
                new CountReorderingMapper(sampleIndex, new NamingConcatFeatureMapper<>(

                        /** Map the from sequence for genotypeIndex=0: */
                        new GenomicContextMapper(indelMappedLength,
                                record -> {
                                    final String fromSequence = record.getSamples(0).getCounts(0).getFromSequence();
                                    return fromSequence.substring(0, Math.min(indelMappedLength, fromSequence.length()));
                                }, true /* no warning if index outside of context, needed since indels have variable lengths */),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(matchesRefMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(originalGobyCountIndexMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(firstBaseMappers),
                        new InverseNormalizationMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(countMappers)),
                        new CeilingNormalizationMapper<>(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(countMappers), 30),
                        new InverseNormalizationMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readIndexMappers)),
                        new GenomicContextMapper(sbiProperties),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(targetAlignedLengthMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(queryAlignedLengthMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(queryPositions),
                        /* NumVariationsInReads for counts not in the best 3: */
                        new DensityMapper("numVariationsInRead",
                                10, sbiProperties,
                                record -> TraversalHelper.forAllSampleCounts(record,
                                        CountInfoOrBuilder::getNumVariationsInReadsList)),

                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(distancesToReadVariations),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(distancesFromStartOfRead),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(distancesFromEndOfRead),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(numVariationsInReadMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readMappingQualityMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(baseQualityMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(bamFlagMappers)
                ));

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
