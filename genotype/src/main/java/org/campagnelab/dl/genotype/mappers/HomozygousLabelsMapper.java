package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * A label that informs if the site is heterozygote (label index 10), or indicates the homozygous genotype called (indices 0-9).
 * Created by rct66 on 12/6/16.
 */
public class HomozygousLabelsMapper extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    private final boolean sortCounts;
    public static final int IS_HETEROZYGOUS_INDEX = 10;
    public static final int NUM_LABELS = IS_HETEROZYGOUS_INDEX + 1;
    private final float epsilon;

    @Override
    public int numberOfLabels() {
        return NUM_LABELS;
    }

    public HomozygousLabelsMapper(boolean sortCounts) {
        this(sortCounts,0);
    }
    public HomozygousLabelsMapper(boolean sortCounts, float epsilon){
            this.sortCounts = sortCounts;
            this.epsilon = epsilon;
        }

        int[] indices = new int[]{0, 0};

        @Override
        public void mapLabels (BaseInformationRecords.BaseInformation record, INDArray labels,int indexOfRecord){
            indices[0] = indexOfRecord;

            for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
                indices[1] = labelIndex;
                labels.putScalar(indices, produceLabel(sortedCountRecord, labelIndex));
            }
        }

        @Override
        public void prepareToNormalize (BaseInformationRecords.BaseInformation record,int indexOfRecord){
            if (sortCounts) {
                sortedCountRecord = sortHelper.sort(record);
            } else {
                sortedCountRecord = record;
            }
        }

        private BaseInformationRecords.BaseInformation sortedCountRecord;
        private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

        @Override
        public float produceLabel (BaseInformationRecords.BaseInformation record,int labelIndex){
            float v = epsilon / (numberOfLabels() + 1);
            record = sortedCountRecord;
            if (!getHomozygous(record)) {
                // this site is heterozygous. The Allele will only be encoded in the other outputs.

                return (labelIndex == IS_HETEROZYGOUS_INDEX/**last index */) ? 1-epsilon : v;
            } else {
                // the site is homozygous.
                if (labelIndex >= record.getSamples(0).getCountsCount()) {
                    // the labelIndex is outside the range of counts in the protobuff.
                    return v;
                } else {
                    // The allele is encoded here:
                    return record.getSamples(0).getCounts(labelIndex).getIsCalled() ? 1-epsilon : v;
                }
            }
        }


    private boolean getHomozygous(BaseInformationRecords.BaseInformation record) {
        record = sortedCountRecord;
        //count number of called alleles
        int numAlleles = 0;
        for (int i = 0; i < record.getSamples(0).getCountsCount(); i++) {
            if (record.getSamples(0).getCounts(i).getIsCalled()) {
                numAlleles++;
            }
        }
        return (numAlleles == 1);
    }


    @Override
    public MappedDimensions dimensions() {
        return new MappedDimensions(numberOfLabels());
    }


}
