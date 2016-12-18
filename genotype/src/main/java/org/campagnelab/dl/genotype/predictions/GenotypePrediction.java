package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;

/**
 * Describes a genotype prediction. Helper method set aggregates individual model predictions.
 * Created by fac2003 on 12/18/16.
 */
public class GenotypePrediction {
    /**
     * Genotype called by the model.
     */
    public String calledGenotype;
    public double overallProbability;
    /**
     * True genotype (if available). Useful for performance evaluation.
     */
    public String trueGenotype;
    /**
     * Probability that a genotype is present. Genotypes are indexed using the goby conventions.
     */
    public double probabilityGenotypeCalled[];
    /**
     * Probability that a genotype is not present. Genotypes are indexed using the goby conventions.
     */
    private double probabilityGenotypeNotCalled[];
    private HomozygousPrediction homozygousPrediction;
    private SingleGenotypePrediction[] singleGenotypePredictions;
    private int numAlleles;

    public void set(HomozygousPrediction homozygousPrediction, SingleGenotypePrediction[] singleGenotypePredictions) {
        this.homozygousPrediction = homozygousPrediction;
        this.singleGenotypePredictions = singleGenotypePredictions;
        if (homozygousPrediction.isHomozygous) {
            calledGenotype = homozygousPrediction.predictedHomozygousGenotype;
            overallProbability = homozygousPrediction.probability;
        } else {
            int genotypeIndex = 0;
            double predProbability = 0;
            StringBuffer hetGenotype = new StringBuffer();
            for (SingleGenotypePrediction singleGenotypePrediction : singleGenotypePredictions) {
                probabilityGenotypeCalled[genotypeIndex] = singleGenotypePrediction.probabilityIsCalled;
                probabilityGenotypeNotCalled[genotypeIndex] = 1 - singleGenotypePrediction.probabilityIsCalled;

                if (singleGenotypePrediction.probabilityIsCalled >= 0.5) {
                    predProbability = singleGenotypePrediction.probabilityIsCalled;
                    numAlleles++;
                    if (hetGenotype.length() > 0) {
                        hetGenotype.append("/");
                    }
                    hetGenotype.append(singleGenotypePrediction.predictedSingleGenotype);
                }
            }
            overallProbability = predProbability / (double) numAlleles;
            calledGenotype = hetGenotype.toString();
        }
    }
}
