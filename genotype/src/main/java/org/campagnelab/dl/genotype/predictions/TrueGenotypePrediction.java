package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;

/**
 * Created by joshuacohen on 2/10/17.
 */
public class TrueGenotypePrediction extends GenotypePrediction {
    public TrueGenotypePrediction(BaseInformationRecords.BaseInformation record, List<Prediction> predictions,
                                  boolean withDistinctAllele, boolean withCombinedLayer, boolean withIsVariant,
                                  double decisionThreshold) {
        GenotypePrediction aggregatePrediction;
        TrueGenotypeOutputLayerPrediction trueGenotypeOutputLayerPrediction;
        MetadataPrediction metadata;
        if (withDistinctAllele) {
            if (withIsVariant) {
                aggregatePrediction = new NumDistinctAlleleWithIsVariantGenotypePrediction(record, decisionThreshold, predictions);
                trueGenotypeOutputLayerPrediction = (TrueGenotypeOutputLayerPrediction) predictions.get(13);
            } else {
                aggregatePrediction = new NumDistinctAlleleGenotypePrediction(record, decisionThreshold, predictions);
                trueGenotypeOutputLayerPrediction = (TrueGenotypeOutputLayerPrediction) predictions.get(12);
            }
            metadata = (MetadataPrediction) predictions.get(11);
        } else if (withCombinedLayer) {
            if (withIsVariant) {
                aggregatePrediction = new CombinedWithIsVariantGenotypePrediction(predictions);
                trueGenotypeOutputLayerPrediction = (TrueGenotypeOutputLayerPrediction) predictions.get(3);
            } else {
                aggregatePrediction = new CombinedGenotypePrediction(predictions);
                trueGenotypeOutputLayerPrediction = (TrueGenotypeOutputLayerPrediction) predictions.get(2);
            }
            metadata = (MetadataPrediction) predictions.get(1);
        } else {
            throw new IllegalArgumentException("The type of aggregate prediction is not recognized.");
        }
//        boolean useAggregate = aggregatePrediction.overallProbability > trueGenotypeOutputLayerPrediction.overallProbability;
        boolean useAggregate = false;
        if (useAggregate) {
            this.overallProbability = aggregatePrediction.overallProbability;
            this.trueGenotype = aggregatePrediction.trueGenotype;
            this.predictedGenotype = aggregatePrediction.predictedGenotype;
            this.isPredictedIndel = aggregatePrediction.isPredictedIndel;
            this.isVariantProbability = aggregatePrediction.isVariantProbability;
        } else {
            this.overallProbability = trueGenotypeOutputLayerPrediction.overallProbability;
            this.trueGenotype = trueGenotypeOutputLayerPrediction.trueGenotype;
            this.predictedGenotype = trueGenotypeOutputLayerPrediction.predictedGenotype;
            this.isPredictedIndel = trueGenotypeOutputLayerPrediction.isPredictedIndel;
            this.isVariantProbability = trueGenotypeOutputLayerPrediction.overallProbability;
        }
        this.isVariant = metadata.isVariant;
        this.isIndel = metadata.isIndel;
        this.index = metadata.index;
        this.referenceGobyIndex = metadata.referenceGobyIndex;
    }
}
