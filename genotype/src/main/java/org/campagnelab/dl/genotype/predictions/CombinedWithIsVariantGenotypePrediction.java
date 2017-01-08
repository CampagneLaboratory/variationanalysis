package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

import java.util.List;

/**
 * Combined with isVariant.
 * Created by fac2003 on 12/23/16.
 */
public class CombinedWithIsVariantGenotypePrediction extends CombinedGenotypePrediction {
    public IsVariantPrediction isVariantPrediction;
    public CombinedWithIsVariantGenotypePrediction(List<Prediction> individualOutputPredictions) {
        super(individualOutputPredictions);
        assert individualOutputPredictions.size()==3:"3 individual predictions are expected";
        this.isVariantPrediction= (IsVariantPrediction) individualOutputPredictions.get(2);
        this.isVariantProbability=isVariantPrediction.probability;
        this.overallProbability=(combinedOutputPrediction.overallProbability+isVariantPrediction.probability)/2;
        this.isVariantProbability=overallProbability;
    }
}
