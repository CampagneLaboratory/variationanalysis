package org.campagnelab.dl.framework.domains.prediction;

/**
 *  Represent a prediction made on a record with a model output. Sub-classes of this class
 *  store interpreted predictions (human understandable, not numeric as produced by neural nets).
 */
public class Prediction {
    /** Index of the record in the source iterator. */
    public int index;
    /**
     * The index of the model output that produced the prediction.
     */
    public int outputIndex;


}
