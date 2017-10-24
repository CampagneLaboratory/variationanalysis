package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;

public class SegmentStatsAccumulator {
    protected ObjectList<String> metrics = new ObjectArrayList();
    protected DoubleArrayList estimates = new DoubleArrayList();
    private SegmentStatsAccumulator accumulators[];

    SegmentStatsAccumulator(SegmentStatsAccumulator... accumulators) {
        this.accumulators = accumulators;
        for (SegmentStatsAccumulator acc : accumulators) {
           metrics.addAll(acc.metricNames());
        }
    }

    public ObjectList<String> metricNames() {

        return metrics;
    }

    public void initializeStats() {

        for (SegmentStatsAccumulator acc : accumulators) {
            acc.initializeStats();
        }
    }


    public void observe(SegmentPrediction fullPred) {
        for (SegmentStatsAccumulator acc : accumulators) {
            acc.observe(fullPred);
        }
    }

    DoubleList estimates() {
        estimates.clear();
        for (SegmentStatsAccumulator acc : accumulators) {
            estimates.addAll(acc.estimates());
        }
        return estimates;
    }
}
