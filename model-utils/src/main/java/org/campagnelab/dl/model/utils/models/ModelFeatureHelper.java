
package org.campagnelab.dl.model.utils.models;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Consumer;

/**
 * Helper class to facilitate predicting output on a network or computation graph.
 * Created by fac2003 on 11/10/16.
 */

public class ModelFeatureHelper<RecordType> {

    private INDArray[] featuresArray;

    public void prepareForNextRecord(Model model, RecordType record,
                                     Consumer<INDArray[]> featuresConsumer,
                                     FeatureMapper... featureMappers) {
        if (featuresArray == null) {
            featuresArray = new INDArray[featureMappers.length];
        }
        if (model instanceof MultiLayerNetwork) {
            MultiLayerNetwork network = (MultiLayerNetwork) model;
            INDArray testFeatures = Nd4j.zeros(1, featureMappers[0].numberOfFeatures());
            featureMappers[0].prepareToNormalize(record, 0);
            featureMappers[0].mapFeatures(record, testFeatures, 0);
            featuresArray[0] = testFeatures;
            featuresConsumer.accept(featuresArray);

        } else if (model instanceof ComputationGraph) {
            ComputationGraph graph = (ComputationGraph) model;
            // INDArray[] testFeatures = new INDArray[featureMappers.length];
            for (int i = 0; i < featureMappers.length; i++) {
                featuresArray[i] = Nd4j.zeros(1, featureMappers[i].numberOfFeatures());
                featureMappers[i].prepareToNormalize(record, 0);
                featureMappers[i].mapFeatures(record, featuresArray[i], 0);
            }

            featuresConsumer.accept(featuresArray);
        } else {
            throw new IllegalArgumentException("model is not of supported type: " + model.getClass().getCanonicalName());
        }
    }

}
