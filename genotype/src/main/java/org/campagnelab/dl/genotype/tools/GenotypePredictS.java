package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.somatic.learning.domains.predictions.SomaticFrequencyPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/12/16.
 */
public class GenotypePredictS extends Predict<BaseInformationRecords.BaseInformation> {


    public static void main(String[] args) {

        Predict predict = new GenotypePredictS();
        predict.parseArguments(args, "PredictS", predict.createArguments());
        predict.execute();
    }

    private AreaUnderTheROCCurve aucLossCalculator;

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        resutsWriter.append("index\ttrueGenotype\tprobabilityA\tprobabilityNo\tcorrectness").append("\n");

    }

    @Override
    protected void initializeStats(String prefix) {
        aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
    }

    //TODO: Implement these
    @Override
    protected void writeOutputStatistics(String prefix, PrintWriter outputWriter) {
        outputWriter.print(aucLossCalculator.evaluateStatistic());
    }

    @Override
    protected void writeOutputHeader(PrintWriter outputWriter) {
        outputWriter.append("auc");
    }

    @Override
    protected void reportStatistics(String prefix) {
        System.out.println("AUC on " + prefix + "=" + aucLossCalculator.evaluateStatistic());
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, List<Prediction> predictionList) {
        // List contains at least one prediction: isSomaticMutation. It may also contain the prediction of
        // somaticFrequency. In the second element, when the model is a computational graph with two outputs.
        BinaryClassPrediction isSomaticMutation = (BinaryClassPrediction) predictionList.get(0);
        String somaticFrequencyText = "";
        if (predictionList.size() >= 2) {
            SomaticFrequencyPrediction somaticFrequency = (SomaticFrequencyPrediction) predictionList.get(1);
            double rmse = somaticFrequency.trueValue == null ? -1 : Math.sqrt(Math.pow(somaticFrequency.trueValue - somaticFrequency.predictedValue, 2));
            somaticFrequencyText += String.format("\t%f\t%f\t%f", somaticFrequency.trueValue, somaticFrequency.predictedValue,
                    rmse);
        }
        String correctness = (isSomaticMutation.predictedLabelYes > isSomaticMutation.predictedLabelNo && isSomaticMutation.trueLabelYes == 1f ||
                isSomaticMutation.predictedLabelNo > isSomaticMutation.predictedLabelYes && isSomaticMutation.trueLabelYes == 0f) ? "correct" : "wrong";

        if (doOuptut(correctness, args(), Math.max(isSomaticMutation.predictedLabelNo, isSomaticMutation.predictedLabelYes))) {
            resultWriter.printf("%d\t%f\t%f\t%f\t%s%s%n", isSomaticMutation.index, isSomaticMutation.trueLabelYes, isSomaticMutation.predictedLabelYes,
                    isSomaticMutation.predictedLabelNo, correctness, somaticFrequencyText);
            if (args().filterAucObservations) {
                aucLossCalculator.observe(isSomaticMutation.predictedLabelYes, isSomaticMutation.trueLabelYes - 0.5);
            }
        }
        //convert true label to the convention used by auc calculator: negative true label=labelNo.
        if (!args().filterAucObservations) {
            aucLossCalculator.observe(isSomaticMutation.predictedLabelYes, isSomaticMutation.trueLabelYes - 0.5);
        }
    }

    /**
     * Apply filters and decide if a prediction should be written to the output.
     *
     * @param correctness
     * @param args
     * @param pMax
     * @return
     */
    protected boolean doOuptut(String correctness, PredictArguments args, double pMax) {
        if (args.correctnessFilter != null) {
            if (!correctness.equals(args.correctnessFilter)) {
                return false;
            }
        }
        if (pMax < args().pFilterMinimum || pMax > args().pFilterMaximum) {
            return false;
        }
        return true;
    }


}