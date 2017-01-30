package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by rct66 on 12/21/16.
 */
public class StatsAccumulatorTest {

    StatsAccumulator acc = new StatsAccumulator();


    @Test
    public void createOutputStatistics() throws Exception {
        observe();
        double[] stats = acc.createOutputStatistics();

        assertEquals("TP is wrong", 1, acc.numTruePositive);
        assertEquals("TN is wrong", 1, acc.numTrueNegative);
        assertEquals("FP is wrong", 1, acc.numFalsePositive);
        assertEquals("FN is wrong", 1, acc.numFalseNegative);


        double[] actual = new double[]{0.5, 0.5, 0.5, 0.5, 2,0.5, Double.NaN, 0.5};
        for (int i = 0; i < stats.length; i++) {
            assertEquals(actual[i], stats[i], 0.00001);
        }
    }


    public void observe() throws Exception {
        acc.initializeStats();
        int nVariants = 0;
        //true negative
        GenotypePrediction pred1 = new GenotypePrediction();
        pred1.trueGenotype = "A/A";
        pred1.predictedGenotype = "A/A";
        pred1.isVariant = false;
        acc.observe(pred1);

        //true positive
        GenotypePrediction pred2 = new GenotypePrediction();
        pred2.trueGenotype = "A/G";
        pred2.predictedGenotype = "G/A";
        pred2.isVariant = true;
        acc.observe(pred2);
        nVariants += 1;

        //false positive
        GenotypePrediction pred3 = new GenotypePrediction();
        pred3.trueGenotype = "T/T";
        pred3.predictedGenotype = "T/C";
        pred3.isVariant = false;
        acc.observe(pred3);

        //false negative
        GenotypePrediction pred4 = new GenotypePrediction();
        pred4.trueGenotype = "A/G";
        pred4.predictedGenotype = "G/G";
        pred4.isVariant = true;
        nVariants += 1;
        acc.observe(pred4);


            assertEquals(1f / nVariants, acc.createOutputStatistics("Recall")[0], 0.1);

    }


}