package org.campagnelab.dl.framework.domains.prediction;

/**
 * Created by joshuacohen on 11/23/16.
 */
public class TimeSeriesPrediction<PredictionType> extends Prediction {
    public PredictionType[] trueLabels;
    public PredictionType[] predictedLabels;
}
