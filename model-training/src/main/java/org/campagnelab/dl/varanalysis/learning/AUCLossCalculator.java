package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.deeplearning4j.earlystopping.scorecalc.ScoreCalculator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.dataset.DataSet;

/**
 * Adapted from DL4J to estimate AUC on the validation set.
 */
public class AUCLossCalculator implements ScoreCalculator<MultiLayerNetwork> {

    private static final int POSITIVE_PROBABILITY_INDEX = 0;
    private static final int NEGATIVE_PROBABILITY_INDEX = 1;
    private DataSetIterator dataSetIterator;
    private boolean average;

    /**
     * Calculate the score (loss function value) on a given data set (usually a test set)
     *
     * @param dataSetIterator Data set to calculate the score for
     * @param average         Whether to return the average (sum of loss / N) or just (sum of loss)
     */
    public AUCLossCalculator(DataSetIterator dataSetIterator, boolean average) {
        this.dataSetIterator = dataSetIterator;
        this.average = average;
        aucCalculator= new AreaUnderTheROCCurve();
    }

    AreaUnderTheROCCurve aucCalculator;

    @Override
    public double calculateScore(MultiLayerNetwork network) {
        dataSetIterator.reset();
        aucCalculator.reset();
        double lossSum = 0.0;
        int exCount = 0;
        while (dataSetIterator.hasNext()) {
            DataSet dataSet = dataSetIterator.next();
            if (dataSet == null) break;
            int nEx = dataSet.getFeatureMatrix().size(0);
            INDArray testPredicted = network.output(dataSet.getFeatures(), false);
            float[] probabilities = testPredicted.getRow(POSITIVE_PROBABILITY_INDEX).data().asFloat();
            float[] labels = dataSet.getLabels().getRow(POSITIVE_PROBABILITY_INDEX).data().asFloat();
            int index = 0;
            for (float prob : probabilities) {

                aucCalculator.observe(prob, labels[index]);
                index+=1;
            }
        }
         return aucCalculator.evaluateStatistic();
    }

    @Override
    public String toString() {
        return "AUCLossCalculator(" + dataSetIterator + ",average=" + average + ")";
    }
}
