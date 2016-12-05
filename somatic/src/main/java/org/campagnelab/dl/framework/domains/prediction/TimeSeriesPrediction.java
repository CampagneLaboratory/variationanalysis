package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

/**
 * Created by joshuacohen on 11/23/16.
 */
public class TimeSeriesPrediction extends Prediction {
    private int[] trueLabels;
    private int[] predictedLabels;

    public TimeSeriesPrediction(int[] trueLabels, int[] predictedLabels) {
        assert trueLabels.length == predictedLabels.length : "Labels should have same length";
        this.trueLabels = trueLabels;
        this.predictedLabels = predictedLabels;
    }

    public TimeSeriesPrediction(int[] trueLabels, INDArray predictedLabels) {
        assert predictedLabels.shape().length == 2 : "Predicted labels should be a 2D array (i.e., just for one time series)";
        assert predictedLabels.shape()[1] == trueLabels.length : "Labels should have same length";
        this.trueLabels = trueLabels;
        this.predictedLabels = Nd4j.argMax(predictedLabels, 1).data().asInt();
    }

    public TimeSeriesPrediction(INDArray trueLabels, INDArray predictedLabels) {
        assert trueLabels.shape().length == 2 : "True labels should be a 2D array (i.e., just for one time series)";
        assert predictedLabels.shape().length == 2 : "Predicted labels should be a 2D array (i.e., just for one time series)";
        assert predictedLabels.shape()[1] == trueLabels.shape()[1] : "Labels should have same length";
        this.trueLabels = getIntArgMaxArray(trueLabels);
        this.predictedLabels = getIntArgMaxArray(predictedLabels);
    }

    public int[] trueLabels() {
        return trueLabels;
    }

    public int[] predictedLabels() {
        return predictedLabels;
    }

    private int[] getIntArgMaxArray(INDArray array) {
        return Nd4j.argMax(array, 1).data().asInt();
    }
}
