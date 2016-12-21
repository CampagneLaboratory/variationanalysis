package org.campagnelab.dl.genotype.tools;


import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.learning.domains.predictions.CombinedPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

/**
 * Prediction class for CombinedGenotype
 *
 * @author Remi Torracinta
 *         Created by rct66 on 12/7/16.
 */
public class PredictGCombined extends PredictG {


    public static void main(String[] args) {

        Predict predict = new PredictGCombined();
        predict.parseArguments(args, "PredictG", predict.createArguments());
        predict.execute();
    }


    @Override
    protected void processPredictions(PrintWriter resultWriter, List<Prediction> predictionList) {
        CombinedPrediction pred = (CombinedPrediction) predictionList.get(0);

        boolean correct = pred.isCorrect();

        String correctness = correct ? "correct" : "wrong";

        if (filterHet((PredictGArguments) args(), pred) && doOuptut(correctness, args(), pred.probability)) {
            resultWriter.printf("%d\t%d\t%s\t%s\t%f\t%s\n",
                    pred.index, (correct ? 1 : 0), pred.trueGenotype, pred.predictedGenotype, pred.probability, correctness);
            if (args().filterMetricObservations) {
                stats.observe(pred, pred.isVariant);
            }
        }
        if (!args().filterMetricObservations) {
            stats.observe(pred, pred.isVariant);
        }

    }


}
