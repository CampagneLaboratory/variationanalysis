package org.campagnelab.dl.framework.performance;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.function.Consumer;
import java.util.function.Predicate;

/**
 * Helper to traverse a record iterator and estimate AUC.
 * Created by fac2003 on 11/3/16.
 */
public class AccuracyHelper {
    public double estimate(DataSetIterator iterator, Model model, int numRecordsForAUC,
                           Consumer<BinaryClassPrediction> doForEachPrediction,
                           Predicate<Integer> stopIfTrue) {

        if (model instanceof ComputationGraph) {
            //This assumes the graph only has one input and predicts the first label:
            final org.deeplearning4j.datasets.iterator.impl.MultiDataSetIteratorAdapter adapter =
                    new org.deeplearning4j.datasets.iterator.impl.MultiDataSetIteratorAdapter(iterator);
            return estimateWithGraph(adapter,
                    (ComputationGraph) model, numRecordsForAUC, doForEachPrediction, stopIfTrue,0);
        }
        throw new RuntimeException("model type not recognized.");
    }



    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, int numRecordsForAUC,
                                    Consumer<BinaryClassPrediction> doForEachPrediction,
                                    Predicate<Integer> stopIfTrue, int outputIndex) {
        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(numRecordsForAUC);
        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;
        BinaryClassPrediction prediction = new BinaryClassPrediction();
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            boolean correct = true;
            for (int predictionIndex = 0; predictionIndex < outputs.length; predictionIndex++) {
                INDArray trueLabels = next.getLabels(outputIndex);
                prediction.trueLabelYes = trueLabels.getDouble(predictionIndex, 1);
                prediction.predictedLabelYes = outputs[outputIndex].getDouble(predictionIndex, 1);
                correct = correct & (prediction.trueLabelYes >= 0.5)?(prediction.predictedLabelYes>=0.5):(prediction.predictedLabelYes<0.5);
                prediction.index = index++;
                doForEachPrediction.accept(prediction);

            }
            if (correct){
                nCorrect++;
            }
            nProcessed++;
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return nCorrect/(double)nProcessed;
    }
}
