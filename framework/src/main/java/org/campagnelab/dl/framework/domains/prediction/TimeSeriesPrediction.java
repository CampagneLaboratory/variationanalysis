package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

/**
 * Created by joshuacohen on 11/23/16.
 */
public class TimeSeriesPrediction extends Prediction {
    public int[] trueLabels;
    public int[] predictedLabels;

    public TimeSeriesPrediction setTrueLabels(int[] trueLabels) {
        if (predictedLabels != null) {
            assert trueLabels.length == predictedLabels.length : "Labels should have same length";
        }
        this.trueLabels = trueLabels;
        return this;
    }

    public TimeSeriesPrediction setTrueLabels(INDArray trueLabels) {
        assert trueLabels.shape().length == 2 : "True labels should be a 2D array (i.e., just for one time series)";
        return setTrueLabels(getIntArgMaxArray(trueLabels));
    }

    public TimeSeriesPrediction setTrueLabels(INDArray allTrueLabels, int labelIdx) {
        assert allTrueLabels.shape().length == 3 : "All true labels should be a 3D array";
        assert labelIdx < allTrueLabels.shape()[0] : "label index is out of bounds";
        return setTrueLabels(allTrueLabels.getRow(labelIdx));
    }

    public TimeSeriesPrediction setPredictedLabels(int[] predictedLabels) {
        if (trueLabels != null) {
            assert trueLabels.length == predictedLabels.length : "Labels should have same length";
        }
        this.predictedLabels = predictedLabels;
        return this;
    }

    public TimeSeriesPrediction setPredictedLabels(INDArray predictedLabels) {
        assert predictedLabels.shape().length == 2 : "Predicted labels should be a 2D array (i.e., just for one time series)";
        return setPredictedLabels(getIntArgMaxArray(predictedLabels));
    }

    public TimeSeriesPrediction setPredictedLabels(INDArray allPredictedLabels, int labelIdx) {
        assert allPredictedLabels.shape().length == 3 : "All predicted labels should be a 3D array";
        assert labelIdx < allPredictedLabels.shape()[0] : "label index is out of bounds";
        return setPredictedLabels(allPredictedLabels.getRow(labelIdx));
    }

    public int[] trueLabels() {
        return trueLabels;
    }

    public int[] predictedLabels() {
        return predictedLabels;
    }

    private int[] getIntArgMaxArray(INDArray array) {
        return Nd4j.argMax(array, 0).data().asInt();
    }
}
