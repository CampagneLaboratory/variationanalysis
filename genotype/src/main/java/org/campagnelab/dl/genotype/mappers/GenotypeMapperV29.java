package org.campagnelab.dl.genotype.mappers;

import org.apache.commons.lang.StringUtils;
import org.campagnelab.dl.framework.mappers.BooleanFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureNameMapper;
import org.campagnelab.dl.framework.mappers.NamedWrapper;
import org.campagnelab.dl.somatic.mappers.*;
import org.campagnelab.dl.somatic.mappers.functional.TraversalHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords.CountInfoOrBuilder;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Properties;

/**
 * V26+ originalGobyCountIndexMapper
 */
public class GenotypeMapperV29 extends GenotypeMapperV11 {


    private FeatureNameMapper<BaseInformationRecords.BaseInformationOrBuilder> delegate;
    //default sampleIndex is zero, adjustable with setter
    private int sampleIndex = 0;
    public GenotypeMapperV29() {
        this(0);
    }
    public GenotypeMapperV29(int sampleIndex) {
        super();
        this.sampleIndex=sampleIndex;
        sortCounts = true;
        withDistinctAlleleCounts = true;
        withCombinedLayer = false;
        MAX_GENOTYPES = 3;
    }
    /**
     * Configure the feature mapper to map a specific sampleIndex
     */
    public void configure(Properties sbiProperties) {


        FeatureNameMapper[] countMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] readIndexMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] readMappingQualityMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] baseQualityMappers = new FeatureNameMapper[MAX_GENOTYPES * 2];
        FeatureNameMapper[] matchesRefMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] isIndelMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] isReadInsertionMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] isReadDeletionMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] firstBaseMappersTo = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] firstBaseMappersFrom = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] numVariationsInReadMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] targetAlignedLengthMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] queryAlignedLengthMappers = new FeatureNameMapper[MAX_GENOTYPES];

        //combined distances for now
        FeatureNameMapper[] distancesToReadVariations = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] bamFlagMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] originalGobyCountIndexMappers = new FeatureNameMapper[MAX_GENOTYPES];
        FeatureNameMapper[] queryPositions = new FeatureNameMapper[MAX_GENOTYPES];




        int genotypeIndex = 0;

        for (int i = 0; i < MAX_GENOTYPES; i++) {
            final int constantGenotypeIndex = genotypeIndex;
            countMappers[i] = (new SingleGenoTypeCountMapper(sampleIndex, i, true));
            readIndexMappers[i] = (new SingleReadIndexCountMapper(sampleIndex, i, true));
            queryPositions[i] = new DensityMapper("queryPosition",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getQueryPositionsList) /*,
                    queryPosition -> (float)(Math.log(queryPosition+1)/Math.log(2))*/);

            matchesRefMappers[i] = (new MatchesReferenceMapper(sampleIndex, i));
            isIndelMappers[i] = new NamedWrapper<BaseInformationRecords.BaseInformationOrBuilder>(
                    new BooleanFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(record ->
                    record.getSamples(sampleIndex).getCounts(constantGenotypeIndex).getIsIndel(), 0)) {
                @Override
                public String getFeatureName(int featureIndex) {
                    return "isIndel_sample="+sampleIndex+"_count="+constantGenotypeIndex;
                }
            };
            isReadInsertionMappers[i] = new NamedWrapper<BaseInformationRecords.BaseInformationOrBuilder>(
                    new BooleanFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(record ->
                            record.getSamples(sampleIndex).getCounts(constantGenotypeIndex).getFromSequence().contains("-"), 0)) {
                @Override
                public String getFeatureName(int featureIndex) {
                    return "isReadInsertion_sample="+sampleIndex+"_count="+constantGenotypeIndex;
                }
            };
            isReadDeletionMappers[i] = new NamedWrapper<BaseInformationRecords.BaseInformationOrBuilder>(
                    new BooleanFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(record ->
                            record.getSamples(sampleIndex).getCounts(constantGenotypeIndex).getToSequence().contains("-"), 0)) {
                @Override
                public String getFeatureName(int featureIndex) {
                    return "isReadDeletion_sample="+sampleIndex+"_count="+constantGenotypeIndex;
                }
            };
            int baseContextLength=10;
            firstBaseMappersTo[i] = new GenomicContextMapper(baseContextLength,

                    record -> {

                        String toSequence = record.getSamples(this.sampleIndex).getCounts(constantGenotypeIndex).getToSequence();
                        int length=Math.min(toSequence.length(),baseContextLength);
                        return StringUtils.rightPad(toSequence.substring(0, length), baseContextLength);

                    });
            firstBaseMappersFrom[i] = new GenomicContextMapper(baseContextLength,

                    record -> {

                        String fromSequence = record.getSamples(this.sampleIndex).getCounts(constantGenotypeIndex).getFromSequence();
                        int length=Math.min(fromSequence.length(),baseContextLength);
                        return StringUtils.rightPad(fromSequence.substring(0, length), baseContextLength);

                    });
            numVariationsInReadMappers[i] = new DensityMapper("numVariationsInRead",
                    10, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotype(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder, BaseInformationRecords.CountInfo::getNumVariationsInReadsList));
            //bin width 1 density mapper that ignores variations outside of caps
            distancesToReadVariations[i] = new DensityMapperCapped("distancesToReadVariations.forward","distancesToReadVariations.reverse",
                    -50,50, sbiProperties,
                    baseInformationOrBuilder ->
                            TraversalHelper.forOneSampleGenotypeBothStrands(sampleIndex, constantGenotypeIndex, baseInformationOrBuilder,
                                    BaseInformationRecords.CountInfo::getDistancesToReadVariationsForwardStrandList,
                                    BaseInformationRecords.CountInfo::getDistancesToReadVariationsReverseStrandList));
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
            bamFlagMappers[i] = new BamFlagMapper(sampleIndex,genotypeIndex);
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
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(originalGobyCountIndexMappers),
                        new NaiveNumAlleleMapper<BaseInformationRecords.BaseInformationOrBuilder>(sampleIndex),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(matchesRefMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(firstBaseMappersTo),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(firstBaseMappersFrom),
                        new InverseNormalizationMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(countMappers)),
                        new InverseNormalizationMapper<BaseInformationRecords.BaseInformationOrBuilder>(
                                new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readIndexMappers)),
                        new GenomicContextMapper(sbiProperties),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(targetAlignedLengthMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(queryAlignedLengthMappers),
                    //    new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(queryPositions),

                        new DensityMapper("numVariationsInRead",
                                10, sbiProperties,
                                record -> TraversalHelper.forAllSampleCounts(record,
                                        CountInfoOrBuilder::getNumVariationsInReadsList)),

                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(distancesToReadVariations),
                          /* NumVariationsInReads for counts not in the best 3: */
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(numVariationsInReadMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(readMappingQualityMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(baseQualityMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(bamFlagMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(isIndelMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(isReadInsertionMappers),
                        new NamingConcatFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>(isReadDeletionMappers)
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
