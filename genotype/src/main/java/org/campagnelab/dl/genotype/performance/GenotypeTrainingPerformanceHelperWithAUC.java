package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.function.Predicate;

/**

 */
public class GenotypeTrainingPerformanceHelperWithAUC extends GenotypeTrainingPerformanceHelper {
    private final GenotypeTrainingPerformanceHelper delegate;
    private double observedScore;
    private double observedAUC;
    AreaUnderTheROCCurve aucCalculator = new AreaUnderTheROCCurve(100000);

    public GenotypeTrainingPerformanceHelperWithAUC(DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor) {
        super(domainDescriptor);
        delegate = new GenotypeTrainingPerformanceHelper(domainDescriptor);
    }

    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, Predicate<Integer> stopIfTrue) {
        final double[] score = new double[1];
        final int[] numMiniBatchesScored = new int[1];

        double delegateReturnValue = delegate.estimateWithGraph(iterator, graph, stopIfTrue, genotypePrediction -> {
                    if (genotypePrediction.isVariant()) {
                        aucCalculator.observe(genotypePrediction.isVariantProbability, genotypePrediction.isCorrect()?1:-1);
                    }
                },
                dsScore -> {
                    score[0] += dsScore;
                    numMiniBatchesScored[0] += 1;
                });

        observedScore = score[0] / (double) numMiniBatchesScored[0];
        observedAUC = aucCalculator.evaluateStatistic();

        return delegateReturnValue;
    }

    public double[] getMetricValues(String... metrics) {
        ObjectArrayList<String> metricsNoScore = ObjectArrayList.wrap(metrics.clone());
        metricsNoScore.remove("score");
        metricsNoScore.remove("AUC");
        metricsNoScore.trim();
        String[] elements = metricsNoScore.toArray(new String[metricsNoScore.size()]);
        DoubleArrayList all = DoubleArrayList.wrap(delegate.getMetricValues(elements));
        for (int i = 0; i < metrics.length; i++) {
            if ("score".equals(metrics[i])) {
                all.add(i, observedScore);
            }
            if ("AUC".equals(metrics[i])) {
                all.add(i, observedAUC);
            }
        }
        return all.toDoubleArray();
    }
}

