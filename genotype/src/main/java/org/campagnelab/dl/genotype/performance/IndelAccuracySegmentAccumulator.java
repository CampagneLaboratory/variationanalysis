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
    private int tp;
    private int fp;
    private int fn;

    @Override
    public ObjectList<String> metricNames() {
        return ObjectArrayList.wrap(new String[]{"indelAccuracy", "#indelCorrect", "#indelPredicted", "#indelTrue",
                "tp", "fp", "fn","indelF1"});
    }

    public void initializeStats() {

        trueIndel = 0;
        predictedIndel = 0;
        trueOrPredictedIndel = 0;
        tp = 0;
        fp = 0;
        fn = 0;
        correct = 0;
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
        tp += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];
                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));

                    return isTrueIndel && trueGenotype.equals(predGenotype);

                });

        fp += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];
                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));
                    // false positives are predicted indels that are not true indels:
                    return !isTrueIndel && predGenotype!=null && predGenotype.contains("-");

                });
        fn += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];
                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));
                    // false negatives are true indels that are not found:
                    return isTrueIndel && !trueGenotype.equals(predGenotype);


                });
        trueOrPredictedIndel += fullPred.numPredictionsWhere(
                (segmentGenotypePrediction, baseIndex) -> {
                    String trueGenotype = segmentGenotypePrediction.trueGenotypes[baseIndex];
                    String predGenotype = segmentGenotypePrediction.predictedGenotypes[baseIndex];

                    boolean isTrueIndel =
                            (trueGenotype != null && trueGenotype.contains("-"));
                    boolean isIndelPredicted =
                            (predGenotype != null && predGenotype.contains("-"));

                    return isTrueIndel || isIndelPredicted;
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
        double f1=2.0*tp/(2.0*tp+fp+fn);
        estimates.clear();
        estimates.add(accuracy);
        estimates.add(correct);
        estimates.add(predictedIndel);
        estimates.add(trueIndel);
        estimates.add(tp);
        estimates.add(fp);
        estimates.add(fn);
        estimates.add(f1);

        return estimates;
    }
}
