package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

import java.util.List;

/**
 * An overall prediction for models that output a combined layer.
 * Homozygous refers only prediction, not the true label.
 */
public class CombinedGenotypePrediction extends GenotypePrediction {

    public CombinedGenotypePrediction(CombinedOutputLayerPrediction outputLayer, MetadataPrediction metaData) {
        this.overallProbability = outputLayer.overallProbability;
        this.trueGenotype = outputLayer.trueGenotype;
        this.predictedGenotype = outputLayer.predictedGenotype;

        this.isVariant = metaData.isVariant;
        this.isIndel = metaData.isIndel;
        this.index = metaData.index;
    }

    public CombinedGenotypePrediction(List<Prediction> individualOutputPredictions) {
        this(/*combined output layer */ (CombinedOutputLayerPrediction)individualOutputPredictions.get(0),
                /*metadata layer */ (MetadataPrediction)individualOutputPredictions.get(1));
    }
}