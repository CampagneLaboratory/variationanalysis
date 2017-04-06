package org.campagnelab.dl.framework.domains.prediction;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.accum.Sum;
import org.nd4j.linalg.api.ops.impl.indexaccum.IMax;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.NDArrayIndex;

/**
 * Created by joshuacohen on 11/23/16.
 */
public class TimeSeriesPrediction extends Prediction {
    public int[] trueLabels;
    public int[] predictedLabels;
    /**
     * Prediction probabilities stored in  probabilities[timeStepIndex][baseIndex]
     */
    public double[][] probabilities;
    Integer predictedSequenceLength;

    public TimeSeriesPrediction() {
        this(null);
    }

    public TimeSeriesPrediction(Integer predictedSequenceLength) {
        this.predictedSequenceLength = predictedSequenceLength;

    }

    public TimeSeriesPrediction setTrueLabels(int[] trueLabels) {
        if (predictedLabels != null) {

            assert trueLabels.length == predictedLabels.length :
                    String.format("Labels should have same length, trueLabels.length=%d predictedLabels.length=%d",
                            trueLabels.length,predictedLabels.length);
        }
        this.trueLabels = trueLabels;
        if (predictedSequenceLength == null) {
            predictedSequenceLength = trueLabels.length;
        }
        return this;
    }

    public double getPredictionProbability(int base, int timeStep) {
        return probabilities[timeStep][base];
    }

    public int numTimeSteps() {

        return probabilities.length;
    }

    public int numBases() {
        return probabilities[0].length;
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
        setProbabilities(predictedLabels);
        return setPredictedLabels(getIntArgMaxArray(predictedLabels));
    }

    public TimeSeriesPrediction setPredictedLabels(INDArray allPredictedLabels, int labelIdx) {
        assert allPredictedLabels.shape().length == 3 : "All predicted labels should be a 3D array";
        assert labelIdx < allPredictedLabels.shape()[0] : "label index is out of bounds";
        int miniBatchSize = allPredictedLabels.size(0);
        setProbabilities(allPredictedLabels, labelIdx);
        return setPredictedLabels(allPredictedLabels.getRow(labelIdx));
    }

    private void setProbabilities(INDArray allPredictedLabels, int labelIdx) {
        int numBases = allPredictedLabels.size(1);
        int numTimeSteps = allPredictedLabels.size(2);
        this.probabilities = new double[numTimeSteps][numBases];
        for (int base = 0; base < numBases; base++) {
            for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) {
                probabilities[timeStep][base] = allPredictedLabels.getDouble(labelIdx, base, timeStep);
            }
        }
    }
    private void setProbabilities(INDArray allPredictedLabels) {
        int numBases = allPredictedLabels.size(0);
        int numTimeSteps = allPredictedLabels.size(1);
        this.probabilities = new double[numTimeSteps][numBases];
        for (int base = 0; base < numBases; base++) {
            for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) {
                probabilities[timeStep][base] = allPredictedLabels.getDouble(base, timeStep);
            }
        }
    }

    public int[] trueLabels() {
        return trueLabels;
    }

    public int[] predictedLabels() {
        return predictedLabels;
    }

    private static int[] getIntArgMaxArray(INDArray array) {
        int maxValidIndex = Nd4j.getExecutioner().exec(new Sum(array), 0).gt(0).sumNumber().intValue();
        INDArray argMax = Nd4j.getExecutioner().exec(new IMax(array), 0);
        return maxValidIndex > 0
                ? argMax.get(NDArrayIndex.all(), NDArrayIndex.interval(0, maxValidIndex)).data().asInt()
                : argMax.data().asInt();
    }
}
