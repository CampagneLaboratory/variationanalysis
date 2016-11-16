package org.campagnelab.dl.varanalysis.learning.domains.predictions;

import org.campagnelab.dl.varanalysis.learning.domains.predictions.Prediction;

/**
 * Represent a binary prediction made on a record.
 */
public class BinaryClassPrediction  extends Prediction {
    /**
     * The true label, 1 when prediction is yes/positive class. 0 otherwise. Null when the true label is unknown.
     */
    public Double trueLabelYes;
    /**
     * Probability that the record is part of the No/negative class.
     */
    public double predictedLabelNo;
    /**
     * Probability that the record is part of the Yes/positive class.
     */
    public double predictedLabelYes;
}