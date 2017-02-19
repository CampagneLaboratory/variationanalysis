package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.mappers.RecordCountSortHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * An interpreter that optionally sorts counts.
 * Created by fac2003 on 12/20/16.
 */
public abstract class SortingCountInterpreter<P extends Prediction> implements PredictionInterpreter<BaseInformationRecords.BaseInformation, P> {
    private boolean sortCounts;
    protected int[] indexPermutation;

    public SortingCountInterpreter(boolean sortCounts) {
        this.sortCounts = sortCounts;
    }

    public BaseInformationRecords.BaseInformation sort(BaseInformationRecords.BaseInformation record) {
        int sortedIndex = 0;
        indexPermutation = new int[record.getSamples(0).getCountsCount()];

        if (sortCounts) {
            sortedCountRecord = sortHelper.sort(record);
            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(0).getCountsList()) {
                indexPermutation[count.getGobyGenotypeIndex()] = sortedIndex++;
            }
        } else {
            sortedCountRecord = record;
            for (BaseInformationRecords.CountInfo count : sortedCountRecord.getSamples(0).getCountsList()) {
                indexPermutation[sortedIndex] = sortedIndex;
                sortedIndex++;
            }
        }
        return sortedCountRecord;
    }

    protected BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

}
