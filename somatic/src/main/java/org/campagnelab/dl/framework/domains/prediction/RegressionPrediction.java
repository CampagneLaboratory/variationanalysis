package org.campagnelab.dl.framework.domains.prediction;

/**
 * Created by fac2003 on 11/12/16.
 */
public class RegressionPrediction extends Prediction {
    public float predictedValue;
    /**
     * The true value, if know, null otherwise.
     */
    public Float trueValue;

}
