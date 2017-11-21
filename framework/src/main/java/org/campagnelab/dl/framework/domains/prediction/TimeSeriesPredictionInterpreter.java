package org.campagnelab.dl.framework.domains.prediction;

import org.campagnelab.goby.util.WarningCounter;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.function.Function;

/**
 * Basic time series prediction interpreter
 * Created by joshuacohen on 11/21/16.
 */
public class TimeSeriesPredictionInterpreter<RecordType> implements PredictionInterpreter<RecordType, TimeSeriesPrediction> {

    static private Logger LOG = LoggerFactory.getLogger(TimeSeriesPredictionInterpreter.class);
WarningCounter warningCounter=new WarningCounter(10);
    private final Function<RecordType, int[]> recordToLabel;
    private final Function<RecordType, Integer> recordToSequenceLength;

    public TimeSeriesPredictionInterpreter(Function<RecordType, int[]> recordToLabel) {
        this(recordToLabel, null);
    }

    public TimeSeriesPredictionInterpreter(Function<RecordType, int[]> recordToLabel,
                                           Function<RecordType, Integer> recordToSequenceLength) {
        this.recordToLabel = recordToLabel;
        this.recordToSequenceLength = recordToSequenceLength;
    }

    @Override
    public TimeSeriesPrediction interpret(RecordType record, INDArray output) {
        int[] trueLabels = recordToLabel.apply(record);
        // true labels have been trimmed to the maximum length supported by the LSTM architecture we are using.
        final Integer predictedSequenceLength = recordToSequenceLength.apply(record);
        if (predictedSequenceLength != trueLabels.length) {
            warningCounter.warn(LOG, "sequence and true labels lengths must agree. " +
                    "Make sure the function recordToSequenceLength accounts for sequence clipping to maxLength.");
        }
        // trim to the maximum length (the LSTM is configured with a max length and we must trim accordingly):
        int maxLength = Math.min(predictedSequenceLength, trueLabels.length);
        TimeSeriesPrediction prediction = new TimeSeriesPrediction(maxLength);

        INDArray trimmedOutput = output.get(NDArrayIndex.all(), NDArrayIndex.interval(0, maxLength));
        prediction.setPredictedLabels(trimmedOutput);
        int[] trimmedTrueLabels = trueLabels.length != maxLength ?
                Arrays.copyOfRange(trueLabels, 0, maxLength) : trueLabels;
        prediction.setTrueLabels(trimmedTrueLabels);
        return prediction;
    }

    @Override
    public TimeSeriesPrediction interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        // note: sequence length ignores masking on the training set.
        TimeSeriesPrediction prediction = new TimeSeriesPrediction();
        prediction.setPredictedLabels(output, exampleIndex);
        prediction.setTrueLabels(trueLabels, exampleIndex);
        return prediction;
    }

    public TimeSeriesPrediction interpret(INDArray trueMasks, INDArray trueLabels, INDArray output, int predictionIndex) {
        assert trueLabels.shape().length == 3 : "True labels should be a 3D array";
        assert trueMasks.shape().length == 2 : "True masks should be a 2D array";
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
