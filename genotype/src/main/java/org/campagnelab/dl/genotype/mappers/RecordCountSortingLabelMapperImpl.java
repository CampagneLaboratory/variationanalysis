package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.somatic.mappers.NoMasksLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 *
 * Mapper that pre-sorts the counts of a record by decreasing allele support. Stores a permutation from the sorted order
 * to  the original count order.
 * Created by fac2003 on 2/21/17.
 */
public abstract class RecordCountSortingLabelMapperImpl extends NoMasksLabelMapper<BaseInformationRecords.BaseInformation> {
    protected int[] indexPermutation;
    protected boolean sortCounts=true;
    protected BaseInformationRecords.BaseInformation sortedCountRecord;
    protected RecordCountSortHelper sortHelper = new RecordCountSortHelper();

    public RecordCountSortingLabelMapperImpl(boolean sortCounts) {
        this.sortCounts = sortCounts;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformation record, int indexOfRecord) {
        int sortedIndex = 0;
        indexPermutation = new int[record.getSamples(this.sampleIndex).getCountsCount()];
        if (sortCounts) {
            sortedCountRecord = sortHelper.sort(0,record);

            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(this.sampleIndex).getCountsList()) {
                indexPermutation[sortedIndex++] = count.getGobyGenotypeIndex();
            }
        } else {
            sortedCountRecord = record;
            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(this.sampleIndex).getCountsList()) {
                indexPermutation[sortedIndex] = sortedIndex;
                sortedIndex++;
            }
        }

    }
}
