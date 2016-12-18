package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * Created by rct66 on 11/12/16.
 */
public class SingleGenotypePrediction extends Prediction {

    /**
     * The allele called/predicted to be present in the sample. e.g., "A" or "A--"
     */
    public String predictedSingleGenotype;
    /**
     * The probability the allele is present, given the model and the data.
     */
    public double probabilityIsCalled;
    /**
     * Whether the allele is present in the ground-truth data.
     */
    public boolean trueIsCalled;


}

