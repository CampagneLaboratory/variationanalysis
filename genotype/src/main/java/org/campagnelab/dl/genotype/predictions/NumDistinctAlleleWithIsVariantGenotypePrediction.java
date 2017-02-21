package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;

/**
 * An interpreter for NumDistinctAllele with a isVariant output.
 * Created by fac2003 on 1/2/17.
 */
public class NumDistinctAlleleWithIsVariantGenotypePrediction extends NumDistinctAlleleGenotypePrediction {
    public IsVariantPrediction isVariantPrediction;
    public NumDistinctAlleleWithIsVariantGenotypePrediction(BaseInformationRecords.BaseInformation record, double decisionThreshold, List<Prediction> predictionList) {
        super(record, decisionThreshold, predictionList);
        assert predictionList.size()==13:"13 predictions are expected";
        this.isVariantPrediction= (IsVariantPrediction) predictionList.get(12);
        this.isVariantProbability=isVariantPrediction.probability;
    }
}
