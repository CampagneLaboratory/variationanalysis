package org.campagnelab.dl.framework.performance;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
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
public class AUCHelper {
    public double estimate(DataSetIterator iterator, Model model, int numRecordsForAUC,
                           Consumer<BinaryClassPrediction> doForEachPrediction,
                           Predicate<Integer> stopIfTrue, PredictionInterpreter interpreter) {
        if (model instanceof MultiLayerNetwork) {
            return estimateWithNet(iterator, (MultiLayerNetwork) model, numRecordsForAUC, doForEachPrediction, stopIfTrue);

        }
        if (model instanceof ComputationGraph) {
            //This assumes the graph only has one input and predicts the first label:
            final org.deeplearning4j.datasets.iterator.impl.MultiDataSetIteratorAdapter adapter =
                    new org.deeplearning4j.datasets.iterator.impl.MultiDataSetIteratorAdapter(iterator);
            return estimateWithGraph(adapter,
                    (ComputationGraph) model, numRecordsForAUC, doForEachPrediction, stopIfTrue, 0,
                    interpreter);
        }
        throw new RuntimeException("model type not recognized.");
    }

    public double estimateWithNet(DataSetIterator iterator, MultiLayerNetwork model, int numRecordsForAUC,
                                  Consumer<BinaryClassPrediction> doForEachPrediction,
                                  Predicate<Integer> stopIfTrue) {
        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(numRecordsForAUC);
        int index = 0;
        int nProcessed = 0;
        BinaryClassPrediction prediction = new BinaryClassPrediction();
        while (iterator.hasNext()) {
            DataSet next = iterator.next();
            INDArray outputs = model.output(next.getFeatures());
            for (int predictionIndex = 0; predictionIndex < next.numExamples(); predictionIndex++) {
                INDArray trueLabels = next.getLabels();
                prediction.trueLabelYes = trueLabels.getDouble(predictionIndex, 1);
                prediction.predictedLabelNo = (float) outputs.getDouble(predictionIndex, 0);
                prediction.predictedLabelYes = (float) outputs.getDouble(predictionIndex, 1);
                aucLossCalculator.observe(prediction.predictedLabelYes, prediction.trueLabelYes - 0.5);
                prediction.index = index++;
                doForEachPrediction.accept(prediction);

            }
            nProcessed += next.numExamples();
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return aucLossCalculator.evaluateStatistic();
    }

    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, int numRecordsForAUC,
                                    Consumer<BinaryClassPrediction> doForEachPrediction,
                                    Predicate<Integer> stopIfTrue, int outputIndex,
                                    PredictionInterpreter interpreter) {
        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(numRecordsForAUC);
        int index = 0;
        int nProcessed = 0;
        BinaryClassPrediction prediction = new BinaryClassPrediction();

        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            int numExamples = next.getFeatures(0).size(0);
            for (int predictionIndex = 0; predictionIndex < numExamples; predictionIndex++) {
                INDArray trueLabels = next.getLabels(outputIndex);
                INDArray output = outputs[outputIndex];

               prediction = (BinaryClassPrediction) interpreter.interpret(trueLabels, output,predictionIndex);
             /*   prediction.trueLabelYes = trueLabels.getDouble(predictionIndex, 0);
                prediction.predictedLabelNo = 1 - outputs[outputIndex].getDouble(predictionIndex, 0);
                prediction.predictedLabelYes = outputs[outputIndex].getDouble(predictionIndex, 0);
*/
                aucLossCalculator.observe(prediction.predictedLabelYes, prediction.trueLabelYes - 0.5);
                prediction.index = index++;
                doForEachPrediction.accept(prediction);

            }
            nProcessed += numExamples;
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return aucLossCalculator.evaluateStatistic();
    }
}
