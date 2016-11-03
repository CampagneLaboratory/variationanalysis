package org.campagnelab.dl.varanalysis.stats;

import htsjdk.samtools.util.Tuple;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;

import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;

/**
 * Helper to traverse a record iterator and estimate AUC.
 * Created by fac2003 on 11/3/16.
 */
public class AUCHelper {

    public double estimate(DataSetIterator iterator,  MultiLayerNetwork model, int numRecordsForAUC,
                           Consumer<Prediction> doForEachPrediction,
                           Predicate<Integer> stopIfTrue) {
        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(numRecordsForAUC);
        int index = 0;
        int nProcessed = 0;
        Prediction prediction = new Prediction();
        while (iterator.hasNext()) {
            DataSet next = iterator.next();
            INDArray outputs = model.output(next.getFeatures());
            for (int predictionIndex = 0; predictionIndex < next.numExamples(); predictionIndex++) {
                INDArray trueLabels = next.getLabels();
                prediction.trueLabelYes = trueLabels.getDouble(predictionIndex, 1);
                prediction.predictedLabelNo = outputs.getDouble(predictionIndex, 0);
                prediction.predictedLabelYes = outputs.getDouble(predictionIndex, 1);
                aucLossCalculator.observe(prediction.predictedLabelYes, prediction.trueLabelYes - 0.5);
                prediction.index=index++;
                doForEachPrediction.accept(prediction);

            }
            nProcessed += next.numExamples();
            if (stopIfTrue.test(nProcessed)) {
                break;
            }

        }
        return aucLossCalculator.evaluateStatistic();
    }
}
