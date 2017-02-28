package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.function.Function;

/**
 * Basic time series prediction interpreter
 * Created by joshuacohen on 11/21/16.
 */
public class TimeSeriesPredictionInterpreter<RecordType> implements PredictionInterpreter<RecordType, TimeSeriesPrediction> {
    private final Function<RecordType, int[]> recordToLabel;

    public TimeSeriesPredictionInterpreter(Function<RecordType, int[]> recordToLabel) {
        this.recordToLabel = recordToLabel;
    }

    @Override
    public TimeSeriesPrediction interpret(RecordType record, INDArray output) {
        int[] trueLabels = recordToLabel.apply(record);
        TimeSeriesPrediction prediction = new TimeSeriesPrediction();
        prediction.setPredictedLabels(output);
        prediction.setTrueLabels(trueLabels);
        return prediction;
    }

    @Override
    public TimeSeriesPrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        TimeSeriesPrediction prediction = new TimeSeriesPrediction();
        prediction.setPredictedLabels(output, predictionIndex);
        prediction.setTrueLabels(trueLabels, predictionIndex);
        return prediction;
    }
}
