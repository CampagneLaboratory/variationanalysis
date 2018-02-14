package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.domains.prediction.BinaryClassPrediction;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedBasePrediction;
import org.campagnelab.dl.somatic.learning.domains.predictions.SomaticFrequencyPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.PrintWriter;
import java.util.List;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/12/16.
 */
public class PredictS extends Predict<BaseInformationRecords.BaseInformation> {


    public static void main(String[] args) {

        Predict predict = new PredictS();
        predict.parseArguments(args, "PredictS", predict.createArguments());
        predict.execute();
    }

    private AreaUnderTheROCCurve aucLossCalculator;

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        boolean hasSomaticFrequency = domainDescriptor.hasOutput("somaticFrequency");
        String somaticFrequencyColumns = hasSomaticFrequency ? "\ttrueSomaticFrequency\tpredictedSomaticFrequency\trmseSomaticFrequency\tpredictedMutatedAllele\ttrueMutatedAllele" : "";
        resutsWriter.append("index\ttrueLabel\tprobabilityYes\tprobabilityNo\tcorrectness" + somaticFrequencyColumns).append("\n");

    }

    @Override
    protected void initializeStats(String prefix) {
        aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
    }

    @Override
    public INDArray getModelOutput(int outputIndex, List<BaseInformationRecords.BaseInformation> records) {
        throw new NotImplementedException();

    }

    private boolean aucCalculated = false;
    private double auc;

    @Override
    protected double[] createOutputStatistics() {
        return new double[]{getAUC(), getAucCiMin(), getAucCiMax()};
    }

    private double getAUC() {
        if (aucCalculated) {
            return auc;
        } else {
            auc = aucLossCalculator.evaluateStatistic();
            aucCalculated = true;
            return auc;
        }
    }

    @Override
    protected String[] createOutputHeader() {
        return new String[]{"auc", "[auc95", "auc95]"};
    }

    @Override
    protected void reportStatistics(String prefix) {
        System.out.println("AUC on " + prefix + "=" + getAUC());
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, BaseInformationRecords.BaseInformation record,
                                      List<Prediction> predictionList) {
        // List contains at least one prediction: isSomaticMutation. It may also contain the prediction of
        // somaticFrequency. In the second element, when the model is a computational graph with two outputs.
        BinaryClassPrediction isSomaticMutation = (BinaryClassPrediction) predictionList.get(0);

        String predictedMutatedAllele = ".";
        String trueMutatedAllele = ".";
        if (isSomaticMutation instanceof IsMutatedBasePrediction) {
            predictedMutatedAllele = ((IsMutatedBasePrediction) isSomaticMutation).predictedMutatedAllele;
            trueMutatedAllele = ((IsMutatedBasePrediction) isSomaticMutation).trueMutatedAllele;
        }

        String somaticFrequencyText = "";
        if (predictionList.size() >= 2) {
            SomaticFrequencyPrediction somaticFrequency = (SomaticFrequencyPrediction) predictionList.get(1);
            double rmse = somaticFrequency.trueValue == null ? -1 : Math.sqrt(Math.pow(somaticFrequency.trueValue - somaticFrequency.predictedValue, 2));
            somaticFrequencyText += String.format("\t%f\t%f\t%f\t%s\t%s", somaticFrequency.trueValue, somaticFrequency.predictedValue,
                    rmse, predictedMutatedAllele, trueMutatedAllele);
        }
        String correctness = (isSomaticMutation.predictedLabelYes > isSomaticMutation.predictedLabelNo && isSomaticMutation.trueLabelYes == 1f ||
                isSomaticMutation.predictedLabelNo > isSomaticMutation.predictedLabelYes && isSomaticMutation.trueLabelYes == 0f) ? "correct" : "wrong";

        if (doOuptut(correctness, args(), Math.max(isSomaticMutation.predictedLabelNo, isSomaticMutation.predictedLabelYes))) {
            resultWriter.printf("%d\t%f\t%f\t%f\t%s%s%n", isSomaticMutation.index, isSomaticMutation.trueLabelYes, isSomaticMutation.predictedLabelYes,
                    isSomaticMutation.predictedLabelNo, correctness, somaticFrequencyText);
            if (args().filterMetricObservations) {
                aucLossCalculator.observe(isSomaticMutation.predictedLabelYes, isSomaticMutation.trueLabelYes - 0.5);
            }
        }
        //convert true label to the convention used by auc calculator: negative true label=labelNo.
        if (!args().filterMetricObservations) {
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


    public double getAucCiMin() {
        return aucLossCalculator.confidenceInterval95()[0];
    }

    public double getAucCiMax() {
        return aucLossCalculator.confidenceInterval95()[1];
    }
}
