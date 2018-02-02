package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.api.Model;
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
    private double observedAUC_F1;

    public GenotypeTrainingPerformanceHelperWithAUC(DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor, Model model) {
        super(domainDescriptor, model);
        delegate = new GenotypeTrainingPerformanceHelper(domainDescriptor, model);
    }

    public double estimateWithGraph(MultiDataSetIterator iterator, ComputationGraph graph, Predicate<Integer> stopIfTrue) {
        final double[] score = new double[1];
        final int[] numMiniBatchesScored = new int[1];

        double delegateReturnValue = delegate.estimateWithGraph(iterator, graph, stopIfTrue, genotypePrediction -> {
                    if (genotypePrediction.isVariant()) {
                        aucCalculator.observe(genotypePrediction.isVariantProbability, genotypePrediction.isCorrect() ? 1 : -1);
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
        metricsNoScore.remove("AUC+F1");
        metricsNoScore.trim();
        String[] elements = metricsNoScore.toArray(new String[metricsNoScore.size()]);
        DoubleArrayList all = DoubleArrayList.wrap(delegate.getMetricValues(elements));

        double F1=0;
//calculate the sum of AUC and F1:
        for (int i = 0; i < metricsNoScore.size(); i++) {
            String metricName = metricsNoScore.get(i);
            if ("F1".equals(metricName)) {
               F1 = all.getDouble(i);
            }
        }
        observedAUC_F1=combine(observedAUC,F1);
        for (int i = 0; i < metrics.length; i++) {
            if ("score".equals(metrics[i])) {
                all.add(i, observedScore);
            }

            if ("AUC".equals(metrics[i])) {
                all.add(i, observedAUC);
            }
            if ("AUC+F1".equals(metrics[i])) {
                all.add(i, observedAUC_F1);
            }
        }
        return all.toDoubleArray();
    }

    private double combine(double auc, double f1) {
        return auc*0.1+f1*0.9;
    }
}

