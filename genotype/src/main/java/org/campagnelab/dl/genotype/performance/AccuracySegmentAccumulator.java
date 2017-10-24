package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;

/**
 * Estimates accuracy of predictions over a segment.
 */
public class AccuracySegmentAccumulator extends SegmentStatsAccumulator {
    int correct = 0;
    int n = 0;

    @Override
    public ObjectList<String> metricNames() {
        return ObjectArrayList.wrap(new String[]{"accuracy","correct","n"});
    }

    public void initializeStats() {
        correct = 0;
        n = 0;
        estimates.clear();
    }

    @Override
    public void observe(SegmentPrediction fullPred) {
        correct += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) ->
                        segmentGenotypePrediction.predictedGenotypes[baseIndex].equals(segmentGenotypePrediction.trueGenotypes[baseIndex]) );
        n += fullPred.numBases();

    }

    @Override
    DoubleList estimates() {
        double accuracy = ((double) correct) / ((double) n);
        estimates.clear();
        estimates.add(accuracy);
        estimates.add(correct);
        estimates.add(n);
        return estimates;
    }
}
