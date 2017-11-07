package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.genotype.predictions.GenotypePrediction;

/**
 * Estimate genotype statistics.
 * Created by rct66 on 12/19/16.
 */
public class StatsAccumulator {


    int numCorrectVariants;
    int numProcessed;
    int numTruePositive;
    int numTrueNegative;
    int numFalsePositive;
    int numFalseNegative;
    int numIndelsCorrect;
    int numSnpsCorrect;
    int numIndelsProcessed;
    int numSnpsProcessed;
    int numIndelsTruePositive;
    int numIndelsFalsePositive;
    int numIndelsFalseNegative;
    int numIndelsTrueNegative;
    int numSnpsTruePositive;
    int numSnpsFalsePositive;
    int numSnpsFalseNegative;
    int numVariants;
    int numIndels;
    int concordantVariants;
    int numVariantsExpected;
    int numTrueOrPredictedVariants;

    int numSnpsTrueNegative;
    int hetCount = 0;
    int homCount = 0;
    int numTrueIndels = 0;
    private int numPredictedIndels = 0;
    private int numIsIndels = 0;
    private int numPredictedSNPs;
    private int numIsSNPs = 0;
    private boolean observedWasFP;
    private boolean observedWasTP;
    private boolean observedWasFN;
    private boolean observedWasTN;

    public void initializeStats() {
        numCorrectVariants = 0;
        numProcessed = 0;
        numTruePositive = 0;
        numTrueNegative = 0;
        numFalsePositive = 0;
        numFalseNegative = 0;
        numVariants = 0;
        numIndels = 0;
        concordantVariants = 0;
        numTrueOrPredictedVariants = 0;
        numIndelsCorrect = 0;
        numSnpsCorrect = 0;
        numIndelsProcessed = 0;
        numSnpsProcessed = 0;
        numSnpsTrueNegative = 0;
        numIndelsTruePositive = 0;
        numIndelsFalsePositive = 0;
        numIndelsFalseNegative = 0;
        numIndelsTrueNegative = 0;
        numSnpsTruePositive = 0;
        numSnpsFalsePositive = 0;
        numSnpsFalseNegative = 0;
        hetCount = 0;
        homCount = 0;
        numPredictedIndels = 0;
        numIsIndels = 0;
    }

    public void observe(GenotypePrediction fullPred) {
        observe(fullPred, fullPred.isVariant(), fullPred.isVariant());
    }

    public void observe(GenotypePrediction fullPred, boolean isTrueVariant, boolean isPredictedVariant) {

        numProcessed++;
        if (isPredictedVariant || isTrueVariant) {
            numTrueOrPredictedVariants += 1;
            concordantVariants += fullPred.isCorrect() ? 1 : 0;
        }
        if (isPredictedVariant) {
            final int size = fullPred.predictedAlleles().size();
            hetCount += (size == 2 ? 1 : 0); //AB
            homCount += (size == 1 ? 1 : 0); //BB
        }
        // estimate FP,TP,FN,TN for SNPs:

        // estimate FP,TP,FN,TN for indels:
        final int foundIndel = isTrueVariant && (fullPred.isIndel()) ? 1 : 0;
        numTrueIndels += foundIndel;
        numIsIndels += fullPred.isIndel() ? 1 : 0;
        numIsSNPs += fullPred.isSnp() ? 1 : 0;
        observedWasFP=false;
        observedWasTP=false;
        observedWasFN=false;
        observedWasTN=false;
        int fp = 0, fn = 0, tp = 0, tn = 0;
        if (fullPred.isCorrect()) {
            if (isTrueVariant) {
                tp = 1;
                observedWasTP=true;
            } else {
                tn = 1;
                observedWasTN=true;
            }
        } else {
            if (isTrueVariant) {
                fn = 1;
                observedWasFN=true;
            } else {
                fp = 1;
                observedWasFP=true;
            }
        }
       if (fullPred.isPredictedIndel() || fullPred.isIndel()) {
            numIndelsTruePositive += tp;
            numIndelsTrueNegative += tn;
            numIndelsFalseNegative += fn;
            numIndelsFalsePositive += fp;
        } else {
            numSnpsTruePositive += tp;
            numSnpsTrueNegative += tn;
            numSnpsFalseNegative += fn;
            numSnpsFalsePositive += fp;
        }
        numPredictedIndels += fullPred.isCorrect() && fullPred.isPredictedIndel() && isTrueVariant ? 1 : 0;
        numVariants += isTrueVariant ? 1 : 0;
        numIndels += fullPred.isIndel() ? 1 : 0;
    }

    public double[] createOutputStatistics() {

        numTrueNegative = numSnpsTrueNegative + numIndelsTrueNegative;
        numTruePositive = numSnpsTruePositive + numIndelsTruePositive;
        numFalseNegative = numSnpsFalseNegative + numIndelsFalseNegative;
        numFalsePositive = numSnpsFalsePositive + numIndelsFalsePositive;

        double recall = ((double)numTruePositive) / ((double) numTruePositive + numFalseNegative);
        double precision = ((double)numTruePositive) / ((double) (numTruePositive + numFalsePositive));

        double F1 = 2 * precision * recall / (precision + recall);
        double indelRecall = ((double)numIndelsTruePositive) / ((double) numIndelsTruePositive + numIndelsFalseNegative);
        double indelPrecision = ((double)numIndelsTruePositive) / ((double) numIndelsTruePositive + numIndelsFalsePositive);
     //        System.out.printf("indels: TP %d FP %d FN %d  TN %d numPredictedIndels: %d %n", numIndelsTruePositive,
     //        numIndelsFalsePositive, numIndelsFalseNegative, numIndelsTrueNegative,numPredictedIndels);
     //        System.out.printf("SNPs:   TP %d FP %d FN %d  TN %d %n", numSnpsTruePositive, numSnpsFalsePositive, numSnpsFalseNegative, numSnpsTrueNegative);
        double indelF1 = 2 * indelPrecision * indelRecall / (indelPrecision + indelRecall);
        double snpRecall = numSnpsTruePositive / ((double) numSnpsTruePositive + numSnpsFalseNegative);
        double snpPrecision = numSnpsTruePositive / ((double) numSnpsTruePositive + numSnpsFalsePositive);
        double snpF1 = 2 * snpPrecision * snpRecall / (snpPrecision + snpRecall);
        double het_hom_ratio = ((double)hetCount)/*AB*/ / (homCount == 0 ? 1 : (double)homCount) /*BB*/;
 //       System.out.printf("het: %d hom: %d ratio: %f %n", hetCount, homCount, het_hom_ratio);
        return new double[]{recall, precision, F1, numVariants,
                indelRecall, indelPrecision, indelF1,
                snpRecall, snpPrecision, snpF1, numIndels,
                het_hom_ratio, numTruePositive, numFalseNegative};
    }

    public double[] createOutputStatistics(String... metrics) {
        double[] estimates = createOutputStatistics();
        double[] values = new double[metrics.length];
        String header[] = createOutputHeader();
        int i = 0;
        for (String metricName : metrics) {
            int metricNameIndex = 0;
            boolean found = false;
            for (metricNameIndex = 0; metricNameIndex < header.length; metricNameIndex++) {
                if (header[metricNameIndex].equals(metricName)) {
                    values[i++] = estimates[metricNameIndex];
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw new RuntimeException("Statistic not found for metric name " + metricName);
            }

        }
        return values;
    }

    public String[] createOutputHeader() {
        return new String[]{"Recall", "Precision", "F1", "NumVariants",
                "Recall_Indels", "Precision_Indels", "F1_Indels",
                "Recall_SNPs", "Precision_SNPs", "F1_SNPs",
                "numIndels", "Het_Hom_Ratio", "TP", "FN"
        };
    }

    public static final int F1_INDEX = 3;

    public void reportStatistics(String prefix) {
        double[] statsArray = createOutputStatistics();
        String[] header = createOutputHeader();
        for (int i = 0; i < Math.min(header.length, statsArray.length); i++) {
            System.out.println(header[i] + "=" + statsArray[i]);
        }
        System.out.printf("Indel TP %d FN %d FP %d TN %d %n", numIndelsTruePositive, numIndelsFalseNegative, numIndelsFalsePositive, numIndelsTrueNegative);
    }

    public void setNumVariantsExpected(int numVariantsExpected) {
        this.numVariantsExpected = numVariantsExpected;
    }

    public boolean observedWasFP() {
        return observedWasFP;
    }
    public boolean observedWasTP() {
        return observedWasTP;
    }
    public boolean observedWasTN() {
        return observedWasTN;
    }
    public boolean observedWasFN() {
        return observedWasFN;
    }
}
