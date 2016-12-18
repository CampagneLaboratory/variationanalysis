package org.campagnelab.dl.framework.performance;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Collections;

/**
 * AUC Calculator. This class was adapted from the BDVAl project. The calculation is straightforward, but has
 * complexity n square, so won't scale with very large number of examples.
 * See http://www.cs.waikato.ac.nz/~remco/roc.pdf for an alternative.
 * Created by fac2003 on 7/15/16.
 *
 * @author Fabien Campagne
 */

public class AreaUnderTheROCCurve {
    static private Logger LOG = LoggerFactory.getLogger(AreaUnderTheROCCurve.class);
    private int maxObservations;
    private boolean clipObservations;
    private DoubleArrayList positiveDecisions;
    private DoubleArrayList negativeDecisions;
    private double estimatedAUC;
    private int numPositive;
    private int numNegative;

    public AreaUnderTheROCCurve() {
        positiveDecisions = new DoubleArrayList();
        negativeDecisions = new DoubleArrayList();
        this.maxObservations = Integer.MAX_VALUE;
        this.clipObservations = false;
    }

    public AreaUnderTheROCCurve(int maxObservations) {
        this();
        this.maxObservations = maxObservations;
        this.clipObservations = true;
    }

    public void reset() {
        positiveDecisions.clear();
        negativeDecisions.clear();
        foundNan = false;
    }

    boolean foundNan = false;

    public void observe(double decisionValue, double label) {
        if (!foundNan && decisionValue != decisionValue) {
            // decision value is NaN:
            LOG.warn("NaN found instead of a decision value. NaN are always interpreted as wrong predictions. ");
            foundNan = true;
        }
        if (label >= 0) {
            positiveDecisions.add(decisionValue);
        } else {
            negativeDecisions.add(decisionValue);
        }
    }

    public double evaluateStatistic() {
        double sum = 0;
        numPositive = 0;
        numNegative = 0;
        if (clipObservations) {
            clipObservations();
        }
        for (final double decisionPositive : positiveDecisions) {
            for (final double decisionNegative : negativeDecisions) {
                sum += decisionPositive > decisionNegative ? 1 : 0;
                sum += decisionPositive == decisionNegative ? 0.5 : 0;
            }
        }

        numPositive = positiveDecisions.size();
        numNegative = negativeDecisions.size();

        final double auc = sum / numPositive / numNegative;
        this.estimatedAUC = auc;
        return auc;
    }

    /**
     * You may call this method after evaluateStatistic() to obtain the 95% confidence interval of the AUC.
     * The confidence interval is estimated using the method of https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Area_Under_an_ROC_Curve.pdf
     *
     * @return The 95% confidence interval.
     */
    public double[] confidenceInterval95() {
        double[] ci = new double[2];
        double z = 1.96f;
        double AUC = estimatedAUC;
        double standardError;
        double Q1 = AUC / (2 - AUC);
        double AUC_square = AUC * AUC;
        double Q2 = 2 * AUC_square / (1 + AUC);
        double N1 = numPositive;
        double N2 = numNegative;

        standardError = Math.sqrt((AUC * (1 - AUC) + (N1 - 1) * (Q1 - AUC_square) + (N2 - 1) * (Q2 - AUC_square))
                / N1 / N2);
        ci[0] = AUC - z * standardError;
        ci[1] = AUC + z * standardError;
        return ci;
    }

    /**
     * Shuffle observations, then clip to the max number. This reduces the precision of the estimate,
     * but ensures computation will finish in a reasonable time (we use a n^2 algorithm).
     */
    private void clipObservations() {
        boolean needToClip = positiveDecisions.size() > maxObservations ||
                negativeDecisions.size() > maxObservations;

        if (needToClip) {
            if (maxObservations < positiveDecisions.size()) {
                Collections.shuffle(positiveDecisions);
                positiveDecisions.size(maxObservations);
            }
            if (maxObservations < negativeDecisions.size()) {
                Collections.shuffle(negativeDecisions);
                negativeDecisions.size(maxObservations);
            }
        }

    }

    public static double evaluateStatistic(final double[] decisionValues, final double[] labels) {
        double sum = 0;
        double numPositive = 0;
        double numNegative = 0;

        final DoubleList truePositiveDecisions = new DoubleArrayList();
        final DoubleList trueNegativeDecisions = new DoubleArrayList();
        for (int i = 0; i < decisionValues.length; i++) {
            if (decisionValues[i] != decisionValues[i]) {
                // decision value is NaN:
                LOG.warn("NaN found instead of a decision value. NaN are always interpreted as wrong predictions. ");
            }
            if (labels[i] >= 0) {
                truePositiveDecisions.add(decisionValues[i]);
            } else {
                trueNegativeDecisions.add(decisionValues[i]);
            }
        }

        for (final double decisionPositive : truePositiveDecisions) {
            for (final double decisionNegative : trueNegativeDecisions) {
                sum += decisionPositive > decisionNegative ? 1 : 0;
                sum += decisionPositive == decisionNegative ? 0.5 : 0;
            }
        }

        numPositive = truePositiveDecisions.size();
        numNegative = trueNegativeDecisions.size();

        final double auc = sum / numPositive / numNegative;
        return auc;
    }
}
