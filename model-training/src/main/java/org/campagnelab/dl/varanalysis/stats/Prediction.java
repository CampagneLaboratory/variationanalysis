package org.campagnelab.dl.varanalysis.stats;

/**
 * Represent a prediction made on a record.
 */
public class Prediction {
    /** Index of the record in the source iterator. */
    public int index;
    public double trueLabelYes;
    public double predictedLabelNo;
    public double predictedLabelYes;
}