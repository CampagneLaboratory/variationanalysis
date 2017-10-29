package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;

/**
 * Estimates accuracy of predictions for indels, over a segment.
 */
public class IndelAccuracySegmentAccumulator extends SegmentStatsAccumulator {
    private int correct;
    private int predictedIndel;
    private int trueIndel;
    private int trueOrPredictedIndel;

    @Override
    public ObjectList<String> metricNames() {
        return ObjectArrayList.wrap(new String[]{"indelAccuracy", "#indelCorrect", "#indelPredicted", "#indelTrue"});
    }

    public void initializeStats() {

        trueIndel=0;
        predictedIndel=0;
        trueOrPredictedIndel=0;
        correct=0;
        estimates.clear();

    }

    @Override
    public void observe(SegmentPrediction fullPred) {
        predictedIndel += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];

                    boolean isIndelPredicted =

                            (predGenotype != null && predGenotype.contains("-"));


                    return isIndelPredicted;
                });

        trueIndel += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];

                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));


                    return isTrueIndel;
                });
        trueOrPredictedIndel += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];

                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));
                    boolean isIndelPredicted =
                            (predGenotype != null && predGenotype.contains("-"));

                    return isTrueIndel||isIndelPredicted;
                });
        correct += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];

                    boolean isIndelCandidate =
                            (trueGenotype != null && trueGenotype.contains("-")) ||
                                    (predGenotype != null && predGenotype.contains("-"));

                    return isIndelCandidate && trueGenotype.equals(predGenotype);
                });

    }

    @Override
    DoubleList estimates() {
        double accuracy = ((double) correct) / ((double) trueOrPredictedIndel);
        estimates.clear();
        estimates.add(accuracy);
        estimates.add(correct);
        estimates.add(predictedIndel);
        estimates.add(trueIndel);
        return estimates;
    }
}
