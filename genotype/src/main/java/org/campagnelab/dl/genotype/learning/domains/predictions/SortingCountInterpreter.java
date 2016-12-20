package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.mappers.RecordCountSortHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * An interprepter that optionally sorts counts.
 * Created by fac2003 on 12/20/16.
 */
public abstract class SortingCountInterpreter<P extends Prediction> implements PredictionInterpreter<BaseInformationRecords.BaseInformation, P> {
    private boolean sortCounts;

    public SortingCountInterpreter(boolean sortCounts) {
        this.sortCounts = sortCounts;
    }

    public BaseInformationRecords.BaseInformation sort(BaseInformationRecords.BaseInformation record) {
        if (sortCounts) {
            sortedCountRecord = sortHelper.sort(record);
        } else {
            sortedCountRecord = record;
        }
        return sortedCountRecord;
    }

    protected BaseInformationRecords.BaseInformation sortedCountRecord;
    private RecordCountSortHelper sortHelper = new RecordCountSortHelper();

}
