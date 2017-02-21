package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Stores meta-data in a virtual label, not used for learning.
 * Created by fac2003 on 12/22/16.
 */
public class MetaDataLabelMapper extends RecordCountSortingLabelMapperImpl {

    public static final int NUM_LABELS = 6;
    /**
     * Stores 1 when the site contains a true variant. Zero otherwise.
     */
    public static final int IS_VARIANT_FEATURE_INDEX = 0;
    /**
     * Stores 1 when the site contains a true indel. Zero otherwise.
     */
    public static final int IS_INDEL_FEATURE_INDEX = 1;
    /**
     * Stores the index of the allele that matches the reference. Index is returned in the original/sbi count index.
     */
    public static final int IS_MATCHING_REF_FEATURE_INDEX = 2;
    /**
     * Permutation from sorted count index to original/sbi count index. For allele with largest count.
     */
    public static final int IS_COUNT1_ORIGINAL_INDEX_FEATURE_INDEX = 3;
    /**
     * Permutation from sorted count index to original/sbi count index. For allele with second largest count.
     */
    public static final int IS_COUNT2_ORIGINAL_INDEX_FEATURE_INDEX = 4;
    /**
     * Permutation from sorted count index to original/sbi count index. For allele with third largest count.
     */
    public static final int IS_COUNT3_ORIGINAL_INDEX_FEATURE_INDEX = 5;

    public MetaDataLabelMapper() {
        super(true);
    }


    @Override
    public int numberOfLabels() {
        return NUM_LABELS;
    }


    int[] indices = new int[]{0, 0};

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformation record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;
        record = sortedCountRecord;
        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformation record, int labelIndex) {
        switch (labelIndex) {
            case IS_VARIANT_FEATURE_INDEX:
                return record.getSamples(0).getIsVariant() ? 1 : 0;
            case IS_INDEL_FEATURE_INDEX:
                final String trueGenotype = record.getTrueGenotype();
                return GenotypeHelper.isIndel(record.getReferenceBase(), trueGenotype) ? 1 : 0;
            case IS_MATCHING_REF_FEATURE_INDEX:
                return calculateReferenceIndex(sortedCountRecord);
            case IS_COUNT1_ORIGINAL_INDEX_FEATURE_INDEX:
                return calculateCountIndex(sortedCountRecord, 0);
            case IS_COUNT2_ORIGINAL_INDEX_FEATURE_INDEX:
                return calculateCountIndex(sortedCountRecord, 1);
            case IS_COUNT3_ORIGINAL_INDEX_FEATURE_INDEX:
                return calculateCountIndex(sortedCountRecord, 2);
            default:
                throw new RuntimeException("No such labelIndex: " + labelIndex);
        }
    }

    private float calculateCountIndex(BaseInformationRecords.BaseInformation record, int sortedCountIndex) {
        return record.getSamples(0).getCounts(sortedCountIndex).getGobyGenotypeIndex();
    }

    /**
     * Determine the original goby count index of the allele matching the reference.
     *
     * @param record
     * @return
     */
    public static int calculateReferenceIndex(BaseInformationRecords.BaseInformation record) {

        BaseInformationRecords.SampleInfo samples = record.getSamples(0);
        for (BaseInformationRecords.CountInfo count : samples.getCountsList()) {
            if (count.getMatchesReference()) {
                return count.getGobyGenotypeIndex();
            }
        }
        return -1;
    }


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }
}

