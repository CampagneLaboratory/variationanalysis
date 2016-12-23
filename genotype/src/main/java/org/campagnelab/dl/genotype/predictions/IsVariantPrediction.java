package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * A prediction that represent the confidence in a correct
 * Created by fac2003 on 12/22/16.
 */
public class IsVariantPrediction extends Prediction {
    /**
     * True when isVariant is predicted for an example.
     */
    public boolean isVariantPredicted;
    /**
     * The probability that isVariant was predicted true.
     */
    public double probability;
    /**
     * The true isVariant from the record.
     */
    public boolean isVariantTruth;
}
