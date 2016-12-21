package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.genotype.mappers.HomozygousLabelsMapper;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.function.Predicate;

/**

 */
public class AccuracyHelper {

    private static final int IS_HOMOZYGOUS_OUTPUT = 0;

    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, Predicate<Integer> stopIfTrue) {


        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;
        BinaryClassPrediction prediction = new BinaryClassPrediction();
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            INDArray[] labels = next.getLabels();

            for (int recordIndex = 0; recordIndex < outputs[IS_HOMOZYGOUS_OUTPUT].rows(); recordIndex++) {
                nProcessed++;
                int homoPredIndex = -1;
                double homoPredMax = -1;

                for (int i = 0; i < outputs[IS_HOMOZYGOUS_OUTPUT].columns(); i++) {
                    double homoPred = outputs[IS_HOMOZYGOUS_OUTPUT].getDouble(recordIndex, i);
                    if (homoPred > homoPredMax) {
                        homoPredIndex = i;
                        homoPredMax = homoPred;
                    }
                }
                if (homoPredIndex != HomozygousLabelsMapper.IS_HETEROZYGOUS_INDEX) {
                    if (labels[IS_HOMOZYGOUS_OUTPUT].getDouble(recordIndex, homoPredIndex) != 0) {
                        nCorrect++;
                    }

                } else {

                    boolean correct = true;
                    for (int predictionIndex = 1; predictionIndex < outputs.length; predictionIndex++) {
                        try {
                            prediction.trueLabelYes = labels[predictionIndex].getDouble(recordIndex, 0);
                        } catch (IndexOutOfBoundsException e) {
                            //System.out.println("me");
                        }
                        prediction.predictedLabelYes = outputs[predictionIndex].getDouble(recordIndex, 0);
                        correct = correct && ((prediction.trueLabelYes >= 0.5) ? (prediction.predictedLabelYes >= 0.5) : (prediction.predictedLabelYes < 0.5));
                        prediction.index = index++;
                    }
                    if (correct) {
                        nCorrect++;
                    }

                }
            }
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return nCorrect/(double)nProcessed;
    }
}
