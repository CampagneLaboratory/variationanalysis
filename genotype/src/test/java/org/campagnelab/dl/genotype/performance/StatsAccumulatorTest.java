package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.genotype.predictions.CombinedGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.junit.Test;

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
        GenotypePrediction pred1 = new GenotypePrediction();
        pred1.trueGenotype = "A/A";
        pred1.predictedGenotype = "A/A";
        pred1.isVariant = false;


        //true positive
        GenotypePrediction pred2 = new GenotypePrediction();
        pred2.trueGenotype = "A/G";
        pred2.predictedGenotype = "G/A";
        pred2.isVariant = true;


        //false positive
        GenotypePrediction pred3 = new GenotypePrediction();
        pred3.trueGenotype = "T/T";
        pred3.predictedGenotype = "T/C";
        pred3.isVariant = false;


        //false negative
        GenotypePrediction pred4 = new GenotypePrediction();
        pred4.trueGenotype = "A/G";
        pred4.predictedGenotype = "G/G";
        pred4.isVariant = true;


        acc.observe(pred1,pred1.isVariant);
        acc.observe(pred2,pred2.isVariant);
        acc.observe(pred3,pred3.isVariant);
        acc.observe(pred4,pred4.isVariant);
assert false;

    }






}