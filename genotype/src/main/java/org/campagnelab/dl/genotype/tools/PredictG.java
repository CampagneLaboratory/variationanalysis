package org.campagnelab.dl.genotype.tools;


import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Remi Torracinta
 *         Created by rct66 on 12/7/16.
 */
public class PredictG extends Predict<BaseInformationRecords.BaseInformation> {

    public static void main(String[] args) {

        Predict predict = new PredictG();
        predict.parseArguments(args, "PredictS", predict.createArguments());
        predict.execute();
    }


    int numCorrect;
    int numProcessed;

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        resutsWriter.append("index\ttrueGenotypeCall\tpredictedGenotypeCall\tcorrectness").append("\n");
    }

    @Override
    protected void initializeStats(String prefix) {
        numCorrect = 0;
        numProcessed = 0;
    }

    //TODO: Implement these
    @Override
    protected void writeOutputStatistics(String prefix, PrintWriter outputWriter) {
        outputWriter.print(numCorrect/(double)numProcessed);
    }

    @Override
    protected void writeOutputHeader(PrintWriter outputWriter) {
        outputWriter.append("accuracy");
    }

    @Override
    protected void reportStatistics(String prefix) {
        System.out.println("Accuracy on " + prefix + "=" + numCorrect/(double)numProcessed);
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, List<Prediction> predictionList) {
        HomozygousPrediction homoPred = (HomozygousPrediction) predictionList.get(0);
        List<Prediction> genoPredList = predictionList.subList(0,11);
        String combinedPredictedGenotype;
        double predProbability;
        combinedPredictedGenotype = homoPred.predictedHomozygousGenotype + ",";
        predProbability = homoPred.probability;
        //now handle the non-homozygous case by concatentating positive single genotype predictions
        if (combinedPredictedGenotype.equals(",")){
            //need running average of single alleles to produce probability
            predProbability = 0;
            int numAlleles = 0;
            StringBuffer nonHomoGenotype = new StringBuffer();
            for (Prediction genoPred : genoPredList) {
                SingleGenotypePrediction singleGenoPred = (SingleGenotypePrediction) genoPred;
                if (singleGenoPred.probability >= 0.5) {
                    predProbability += singleGenoPred.probability;
                    numAlleles++;
                    nonHomoGenotype.append(singleGenoPred.predictedSingleGenotype + ",");
                }
            }
            combinedPredictedGenotype = nonHomoGenotype.toString();
            predProbability = predProbability/(double)numAlleles;
        }
        boolean correct = combinedPredictedGenotype.equals(homoPred.trueGenotype);


        //remove dangling commas
        String correctness = correct ? "correct" : "wrong";
        String trueGenotype = homoPred.trueGenotype.substring(0,homoPred.trueGenotype.length()-1);
        String predictedGenotype = combinedPredictedGenotype.substring(0,combinedPredictedGenotype.length()-1);

        if (doOuptut(correctness, args(), predProbability)) {
            resultWriter.printf("%d\t%s\t%s\t%b", homoPred.index, trueGenotype, predictedGenotype, correctness);
            if (args().filterMetricObservations) {
                numProcessed++;
                if (correct) {
                    numCorrect++;
                };
            }
        }
        //convert true label to the convention used by auc calculator: negative true label=labelNo.
        if (!args().filterMetricObservations) {
            numProcessed++;
            if (correct) {
                numCorrect++;
            };
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
