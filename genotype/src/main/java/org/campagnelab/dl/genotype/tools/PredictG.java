package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousPrediction;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
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

    @Override
    public PredictArguments createArguments() {
        return new PredictGArguments();
    }

    public static void main(String[] args) {

        Predict predict = new PredictG();
        predict.parseArguments(args, "PredictG", predict.createArguments());
        predict.execute();
    }


    int numCorrect;
    int numProcessed;
    int numTruePositive;
    int numTrueNegative;
    int numFalsePositive;
    int numFalseNegative;

    double accuracy;
    double recall;
    double precision;
    double F1;

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        resutsWriter.append("index\tpredictionCorrect01\ttrueGenotypeCall\tpredictedGenotypeCall\tprobabilityIsCalled\tcorrectness").append("\n");
    }

    @Override
    protected void initializeStats(String prefix) {
        numCorrect = 0;
        numProcessed = 0;
        numTruePositive = 0;
        numTrueNegative = 0;
        numFalsePositive = 0;
        numFalseNegative = 0;
    }



    @Override
    protected double[] createOutputStatistics() {
        accuracy = numCorrect/(double)numProcessed;
        recall = numTruePositive/((double)(numTruePositive+numFalseNegative));
        precision = numTruePositive/((double)(numTruePositive+numFalsePositive));
        F1 = precision*recall/(precision+recall);
        return new double[]{accuracy,recall,precision,F1};
    }

    @Override
    protected String[] createOutputHeader() {
        return new String[]{"accuracy","sensitivity(recall)","PPV(precision)","F1"};
    }

    @Override
    protected void reportStatistics(String prefix) {
        System.out.println("Accuracy on " + prefix + "=" + accuracy);
        System.out.println("Recall on " + prefix + "=" + recall);
        System.out.println("Precision on " + prefix + "=" + precision);
        System.out.println("F1 on " + prefix + "=" + F1);
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, List<Prediction> predictionList) {
        HomozygousPrediction homoPred = (HomozygousPrediction) predictionList.get(0);
        List<Prediction> genoPredList = predictionList.subList(1,11);

        GenotypePrediction fullPred = new GenotypePrediction();
        fullPred.set(homoPred,genoPredList.toArray(new SingleGenotypePrediction[genoPredList.size()]));

        boolean correct = fullPred.isCorrect();
        //remove dangling commas
        String correctness = correct ? "correct" : "wrong";

        if (filterHet((PredictGArguments) args(),fullPred) && doOuptut(correctness, args(), fullPred.overallProbability)) {
            resultWriter.printf("%d\t%d\t%s\t%s\t%f\t%s\n",
                    homoPred.index, (correct?1:0), fullPred.trueGenotype, fullPred.calledGenotype, fullPred.overallProbability, correctness);
            if (args().filterMetricObservations) {
                numProcessed++;
                if (correct) {
                    numCorrect++;
                    if (homoPred.isVariant) {
                        numTruePositive++;
                    } else {
                        numTrueNegative++;
                    }
                } else {
                    if (homoPred.isVariant) {
                        numFalseNegative++;
                    } else {
                        numFalsePositive++;
                    }
                }
            }
        }
        //convert true label to the convention used by auc calculator: negative true label=labelNo.
        if (!args().filterMetricObservations) {
            numProcessed++;
            if (correct) {
                numCorrect++;
                if (homoPred.isVariant) {
                    numTruePositive++;
                } else {
                    numTrueNegative++;
                }
            } else {
                if (homoPred.isVariant) {
                    numFalseNegative++;
                } else {
                    numFalsePositive++;
                }
            }
        }
    }

    private boolean filterHet(PredictGArguments args, GenotypePrediction fullPred) {
        Set<String> alleles = fullPred.alleles();
        switch (args.showFilter) {
            case HET:
                return (alleles.size()==2);
            case HOM:
                return alleles.size()==1;
            default:
                return true;
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
