package org.campagnelab.dl.framework.domains.prediction;

/**
 * Represent a binary prediction made on a record.
 */
public class MultiClassPrediction extends Prediction {

    /**
     * Probability of each class. Should sum to 1.
     */
    public int trueClassIndex;
    /**
     * Probability that the record is part of the No/negative class.
     */
    public double[] predictedProbabilities;

}
