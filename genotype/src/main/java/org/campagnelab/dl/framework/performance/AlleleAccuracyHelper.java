package org.campagnelab.dl.framework.performance;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.function.Predicate;

/**
 * Helper to traverse a record iterator and estimate AUC.
 * Created by fac2003 on 11/3/16.
 */
public class AlleleAccuracyHelper {

    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, Predicate<Integer> stopIfTrue) {
        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;
        BinaryClassPrediction prediction = new BinaryClassPrediction();
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            INDArray[] labels = next.getLabels();

            for (int recordIndex = 0; recordIndex < outputs[0].rows(); recordIndex++){
                for (int predictionIndex = 1; predictionIndex < outputs.length; predictionIndex++) {
                    nProcessed++;
                    prediction.trueLabelYes = labels[predictionIndex].getDouble(recordIndex, 0);
                    prediction.predictedLabelYes = outputs[predictionIndex].getDouble(recordIndex, 0);
                    if ((prediction.trueLabelYes >= 0.5)?(prediction.predictedLabelYes>=0.5):(prediction.predictedLabelYes<0.5)){
                        nCorrect++;
                    }
                    prediction.index = index++;
                }
            }
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return nCorrect/(double)nProcessed;
    }
}
