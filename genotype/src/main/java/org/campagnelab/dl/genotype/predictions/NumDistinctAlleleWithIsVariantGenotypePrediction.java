package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

import java.util.List;

/**
 * An interpreter for NumDistinctAllele with a isVariant output.
 * Created by fac2003 on 1/2/17.
 */
public class NumDistinctAlleleWithIsVariantGenotypePrediction extends NumDistinctAlleleGenotypePrediction {
    public IsVariantPrediction isVariantPrediction;
    public NumDistinctAlleleWithIsVariantGenotypePrediction(double decisionThreshold, List<Prediction> predictionList) {
        super(decisionThreshold, predictionList);
        assert predictionList.size()==13:"13 predictions are expected";
        this.isVariantPrediction= (IsVariantPrediction) predictionList.get(12);
        this.isVariantProbability=isVariantPrediction.probability;
    }
}
