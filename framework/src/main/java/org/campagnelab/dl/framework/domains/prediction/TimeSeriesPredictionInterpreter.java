package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.indexing.NDArrayIndex;

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
        // TODO: Make sure trimming the output doesn't break anything
        INDArray trimmedOutput = output.get(NDArrayIndex.all(), NDArrayIndex.interval(0, trueLabels.length));
        prediction.setPredictedLabels(trimmedOutput);
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

    public TimeSeriesPrediction interpret(INDArray trueMasks, INDArray trueLabels, INDArray output, int predictionIndex) {
        assert trueLabels.shape().length == 3 : "True labels should be a 3D array";
        assert trueMasks.shape().length == 3 : "True masks should be a 3D array";
        assert output.shape().length == 3 : "True labels should be a 3D array";
        assert predictionIndex < trueMasks.shape()[0] : "prediction index is out of bounds for true masks";
        assert predictionIndex < trueLabels.shape()[0] : "prediction index is out of bounds for true labels";
        assert predictionIndex < output.shape()[0] : "prediction index is out of bounds for output";
        int maxTrueLabels = trueMasks.getRow(predictionIndex).gt(0).sumNumber().intValue();
        INDArray trimmedTrueLabels = trueLabels.get(NDArrayIndex.all(), NDArrayIndex.all(),
                NDArrayIndex.interval(0, maxTrueLabels));
        INDArray trimmedOutput = output.get(NDArrayIndex.all(), NDArrayIndex.all(),
                NDArrayIndex.interval(0, maxTrueLabels));
        return interpret(trimmedTrueLabels, trimmedOutput, predictionIndex);
    }
}
