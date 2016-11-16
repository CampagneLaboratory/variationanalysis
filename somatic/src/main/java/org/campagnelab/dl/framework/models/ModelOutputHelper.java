
package org.campagnelab.dl.framework.models;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.DataSet;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;
import java.util.Iterator;

/**
 * Helper class to facilitate predicting output on a network or computation graph.
 * Created by fac2003 on 11/10/16.
 */

public class ModelOutputHelper<RecordType> {
    private INDArray[] resultGraph;

    /**
     *
     * @param model
     * @param iterator Must be of type Iterator<DataSet> or Iterator<MultiDataSet>.
     */
    @SuppressWarnings("unchecked")
    public void predictForNext(Model model, Iterator iterator) {
        // we cannot check the casts below because we cannot test for generics with instanceof:
        if (model instanceof MultiLayerNetwork) {
            predictForNext((MultiLayerNetwork) model, (Iterator<DataSet>) iterator);
        } else if (model instanceof ComputationGraph) {
            predictForNext((ComputationGraph) model, (Iterator<MultiDataSet>) iterator);
        }
    }

    public void predictForNextRecord(Model model, RecordType record, FeatureMapper... featureMappers) {

          if (model instanceof MultiLayerNetwork) {
            MultiLayerNetwork network = (MultiLayerNetwork) model;
            INDArray testFeatures = Nd4j.zeros(1, featureMappers[0].numberOfFeatures());
            featureMappers[0].prepareToNormalize(record, 0);
            featureMappers[0].mapFeatures(record, testFeatures, 0);
            Arrays.fill(resultGraph, null);
            resultGraph[0] = network.output(testFeatures, false);
        } else if (model instanceof ComputationGraph) {
            ComputationGraph graph = (ComputationGraph) model;
            INDArray[] testFeatures = new INDArray[featureMappers.length];
            for (int i = 0; i < featureMappers.length; i++) {
                testFeatures[i] = Nd4j.zeros(1, featureMappers[i].numberOfFeatures());
                featureMappers[i].prepareToNormalize(record, 0);
                featureMappers[i].mapFeatures(record, testFeatures[i], 0);
            }
            resultGraph = graph.output(false, testFeatures);
        }else {
              throw new IllegalArgumentException("model is not of supported type: "+model.getClass().getCanonicalName());
          }
    }

    public void predictForNext(ComputationGraph graph, Iterator<MultiDataSet> iterator) {
        resultGraph = graph.output(false, iterator.next().getFeatures());
    }

    public void predictForNext(MultiLayerNetwork network, Iterator<DataSet> iterator) {
        Arrays.fill(resultGraph, null);
        resultGraph[0] = network.output(iterator.next().getFeatures(), false);
    }

    public INDArray getOutput(int outputIndex) {
        return resultGraph[outputIndex];
    }
}
