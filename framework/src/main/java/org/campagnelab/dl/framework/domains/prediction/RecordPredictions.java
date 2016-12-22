package org.campagnelab.dl.framework.domains.prediction;

import java.util.List;

/**
 * Created by fac2003 on 12/22/16.
 */
public class RecordPredictions<RecordType> {
    public RecordPredictions(RecordType record, List<Prediction> predictions) {
        this.record = record;
        this.predictions = predictions;
    }

    public RecordType record;
    public   List<Prediction> predictions;
}
