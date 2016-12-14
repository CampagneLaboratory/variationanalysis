package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

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
        resutsWriter.append("index\ttrueGenotypeCall\tpredictedGenotypeCall\tprobabilityIsCalled\tcorrectness").append("\n");
    }

    @Override
    protected void initializeStats(String prefix) {
        numCorrect = 0;
        numProcessed = 0;
    }



    @Override
    protected double[] createOutputStatistics() {
        return new double[]{numCorrect/(double)numProcessed};
    }

    @Override
    protected String[] createOutputHeader() {
        return new String[]{"accuracy"};
    }

    @Override
    protected void reportStatistics(String prefix) {
        System.out.println("Accuracy on " + prefix + "=" + numCorrect/(double)numProcessed);
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, List<Prediction> predictionList) {
        HomozygousPrediction homoPred = (HomozygousPrediction) predictionList.get(0);
        List<Prediction> genoPredList = predictionList.subList(1,11);

        //use format just for stats output writing, correctness determined with string sets.
        String predictedGenotypeFormat;
        //this set will be compared against true genotype set to check for prediction correctness
        Set<String> predictedGenotype = new ObjectArraySet<>();
        double predProbability;

        //first try homozygous
        predictedGenotype.add(homoPred.predictedHomozygousGenotype);
        predictedGenotypeFormat = homoPred.predictedHomozygousGenotype + "/" + homoPred.predictedHomozygousGenotype ;
        predProbability = homoPred.probability;


        //now handle the non-homozygous case by concatentating positive single genotype predictions
        if (predictedGenotypeFormat.equals("/")){
            //need running average of single alleles to produce probabilityIsCalled
            predProbability = 0;
            int numAlleles = 0;
            StringBuffer nonHomoGenotype = new StringBuffer();
            for (Prediction genoPred : genoPredList) {
                SingleGenotypePrediction singleGenoPred = (SingleGenotypePrediction) genoPred;
                if (singleGenoPred.probabilityIsCalled >= 0.5) {
                    predProbability += singleGenoPred.probabilityIsCalled;
                    numAlleles++;
                    nonHomoGenotype.append(singleGenoPred.predictedSingleGenotype + "/");
                    predictedGenotype.add(singleGenoPred.predictedSingleGenotype);

                }
            }
            predictedGenotypeFormat = nonHomoGenotype.toString();
            try {
                predictedGenotypeFormat = predictedGenotypeFormat.substring(0, predictedGenotypeFormat.length() - 1);
            } catch ( StringIndexOutOfBoundsException e ){
                System.out.println("example with no calls found at index" + homoPred.index);
            }
            predProbability = predProbability/(double)numAlleles;
        }
        boolean correct = predictedGenotype.equals(homoPred.trueGenotype);


        //remove dangling commas
        String correctness = correct ? "correct" : "wrong";

        if (doOuptut(correctness, args(), predProbability)) {
            resultWriter.printf("%d\t%s\t%s\t%f\t%s\n", homoPred.index, homoPred.trueGenotypeFormat, predictedGenotypeFormat, predProbability, correctness);
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
