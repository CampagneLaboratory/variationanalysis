package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.prediction.Prediction;

import java.util.ArrayList;

/**
 * Represent the "genotype" output predictions in human readable form.
 */
public class SegmentGenotypePrediction extends Prediction {
    // predictedGenotypes[siteIndex] contains the predicted genotype (e.g., G/-, T/T, etc.)
    public String[] predictedGenotypes;

    // trueGenotypes[siteIndex] contains the true genotype (e.g., G/-, T/T, etc.), obtained from the segment record.
    public String[] trueGenotypes;

    // probability of prediction, one per base predicted
    public float[] probabilities;
    public int length;

    public int numBases() {
        return length;
        //return predictedGenotypes.length;
    }

    @Override
    public String toString() {
        return String.format("predicted: %s\n" + "true:      %s\n", ObjectArrayList.wrap(predictedGenotypes),
                ObjectArrayList.wrap(trueGenotypes));
    }
}