package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * Created by fac2003 on 7/2/17.
 */
public class SoftmaxGenotypePrediction extends Prediction {
    public int predictedGenotypeIndex;
    public double probability;
    public String trueGenotype;
    public int trueGenotypeIndex;
    public int numBits;
}
