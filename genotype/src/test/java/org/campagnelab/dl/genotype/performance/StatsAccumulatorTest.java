package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.genotype.learning.domains.predictions.CombinedPrediction;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by rct66 on 12/21/16.
 */
public class StatsAccumulatorTest {

    StatsAccumulator acc = new StatsAccumulator();


    @Test
    public void createOutputStatistics() throws Exception {
        observe();
        double[] stats = acc.createOutputStatistics();

        assert acc.numTruePositive == 1;
        assert acc.numTrueNegative == 1;
        assert acc.numFalsePositive == 1;
        assert acc.numFalseNegative == 1;


        double[] actual = new double[]{0.5,0.5,0.5,.25,2};
        for (int i = 0; i < stats.length; i++){
            assert Math.abs(stats[i] - actual[i]) < 0.00001;
        }
    }


    public void observe() throws Exception {
        acc.initializeStats();

        //true negative
        CombinedPrediction pred1 = new CombinedPrediction();
        pred1.trueGenotypeFormat = "A/A";
        pred1.predictedGenotype = "A/A";
        pred1.isVariant = false;


        //true positive
        CombinedPrediction pred2 = new CombinedPrediction();
        pred2.trueGenotypeFormat = "A/G";
        pred2.predictedGenotype = "G/A";
        pred2.isVariant = true;


        //false positive
        CombinedPrediction pred3 = new CombinedPrediction();
        pred3.trueGenotypeFormat = "T/T";
        pred3.predictedGenotype = "T/C";
        pred3.isVariant = false;


        //false negative
        CombinedPrediction pred4 = new CombinedPrediction();
        pred4.trueGenotypeFormat = "A/G";
        pred4.predictedGenotype = "G/G";
        pred4.isVariant = true;


        acc.observe(pred1,pred1.isVariant);
        acc.observe(pred2,pred2.isVariant);
        acc.observe(pred3,pred3.isVariant);
        acc.observe(pred4,pred4.isVariant);


    }






}