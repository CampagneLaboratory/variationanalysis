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
    public void testOutputStatisticsSNPs() throws Exception {
        observe(false);
        double[] stats = acc.createOutputStatistics();

        assertEquals("TP is wrong", 1, acc.numTruePositive);
        assertEquals("TN is wrong", 1, acc.numTrueNegative);
        assertEquals("FP is wrong", 1, acc.numFalsePositive);
        assertEquals("FN is wrong", 1, acc.numFalseNegative);

        assertEquals("TP is wrong", 1, acc.numSnpsTruePositive);
        assertEquals("TN is wrong", 1, acc.numSnpsTrueNegative);
        assertEquals("FP is wrong", 1, acc.numSnpsFalsePositive);
        assertEquals("FN is wrong", 1, acc.numSnpsFalseNegative);
        //  return new double[]{recall, precision, F1, numVariants,
        //        indelRecall, indelPrecision, indelF1,
        //      snpRecall, snpPrecision, snpF1, numIndels,
        //    het_hom_ratio, numTruePositive, numTrueNegative};

        String header[] = acc.createOutputHeader();
        double[] actual = new double[]{0.5, 0.5, 0.5, 2.0,
                Double.NaN, Double.NaN, Double.NaN,
                0.5, 0.5, 0.5,  0,
                1, 1, 1};
        for (int i = 0; i < stats.length; i++) {
            assertEquals("wrong stat at index " + i + " " + header[i], actual[i], stats[i], 0.00001);
        }
    }

    @Test
    public void testOutputStatisticsIndels() throws Exception {
        observe(true);
        double[] stats = acc.createOutputStatistics();

        assertEquals("TP is wrong", 1, acc.numTruePositive);
        assertEquals("TN is wrong", 1, acc.numTrueNegative);
        assertEquals("FP is wrong", 1, acc.numFalsePositive);
        assertEquals("FN is wrong", 1, acc.numFalseNegative);

        assertEquals("TP is wrong", 1, acc.numIndelsTruePositive);
        assertEquals("TN is wrong", 1, acc.numIndelsTrueNegative);
        assertEquals("FP is wrong", 1, acc.numIndelsFalsePositive);
        assertEquals("FN is wrong", 1, acc.numIndelsFalseNegative);

        //  return new double[]{recall, precision, F1, numVariants,
        //        indelRecall, indelPrecision, indelF1,
        //      snpRecall, snpPrecision, snpF1, numIndels,
        //    het_hom_ratio, numTruePositive, numTrueNegative};


    String header[] = acc.createOutputHeader();
    double[] actual = new double[]{0.5, 0.5, 0.5,2.0,
            0.5, 0.5, 0.5,
            Double.NaN, Double.NaN, Double.NaN, 4.0,
            2.0, 1.0, 1.0};
        for(
    int i = 0;
    i<stats.length;i++)

    {
        assertEquals("wrong stat at index " + i + " " + header[i], actual[i], stats[i], 0.00001);
    }

}

    public void observe(boolean isIndel) throws Exception {
        acc.initializeStats();
        int nVariants = 0;
        //true negative
        GenotypePrediction pred1 = new GenotypePrediction();
        pred1.trueGenotype = "A/A" + (isIndel ? "-" : "");
        pred1.predictedGenotype = "A/A" + (isIndel ? "-" : "");
        pred1.isVariant = false;
        pred1.isIndel = isIndel;
        pred1.isPredictedIndel=isIndel;
        acc.observe(pred1);

        //true positive
        GenotypePrediction pred2 = new GenotypePrediction();
        pred2.trueGenotype = "A" + (isIndel ? "-" : "") + "/G";
        pred2.predictedGenotype = "G/A" + (isIndel ? "-" : "");
        pred2.isVariant = true;
        pred2.isIndel = isIndel;
        pred1.isPredictedIndel=isIndel;
        acc.observe(pred2);
        nVariants += 1;

        //false positive
        GenotypePrediction pred3 = new GenotypePrediction();
        pred3.trueGenotype = "T/T" + (isIndel ? "-" : "");
        pred3.predictedGenotype = "T/C" + (isIndel ? "-" : "");
        pred3.isVariant = false;
        pred3.isIndel = isIndel;
        pred1.isPredictedIndel=isIndel;
        acc.observe(pred3);

        //false negative
        GenotypePrediction pred4 = new GenotypePrediction();
        pred4.trueGenotype = "A/G" + (isIndel ? "-" : "");
        pred4.predictedGenotype = "G/G" + (isIndel ? "-" : "");
        pred4.isVariant = true;
        pred4.isIndel = isIndel;
        pred1.isPredictedIndel=isIndel;
        nVariants += 1;
        acc.observe(pred4);


        //      assertEquals(1f / nVariants, acc.createOutputStatistics("Recall")[0], 0.1);

    }


}