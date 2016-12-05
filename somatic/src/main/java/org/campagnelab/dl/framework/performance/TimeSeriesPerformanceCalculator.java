package org.campagnelab.dl.framework.performance;


import org.campagnelab.dl.framework.domains.prediction.TimeSeriesPrediction;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by joshuacohen on 11/28/16.
 */
public class TimeSeriesPerformanceCalculator {
    private Counter truePositives;
    private Counter falsePositives;
    private Counter falseNegatives;
    private Set<Integer> neverPredictedLabels;
    private Set<Integer> neverAppearedLabels;
    private int correctPredictions;
    private int totalPredictions;
    private Map<Integer, Double> precision;
    private Map<Integer, Double> recall;
    private ConfusionMatrix confusionMatrix;
    private List<Integer> allLabels;
    private Map<String, Object> stats;


    public TimeSeriesPerformanceCalculator(List<Integer> allLabels) {
        truePositives = new Counter();
        falsePositives = new Counter();
        falseNegatives = new Counter();
        neverPredictedLabels = new HashSet<>();
        neverAppearedLabels = new HashSet<>();
        precision = new HashMap<>();
        recall = new HashMap<>();
        confusionMatrix = new ConfusionMatrix(allLabels);
        totalPredictions = 0;
        stats = new HashMap<>();
        this.allLabels = allLabels;
    }

    public TimeSeriesPerformanceCalculator(int numLabels) {
        this(IntStream.range(0, numLabels).boxed().collect(Collectors.toList()));
    }

    public void addTimeSeries(TimeSeriesPrediction timeSeries) {
        for (int i = 0; i < timeSeries.trueLabels().length; i++) {
            totalPredictions++;
            confusionMatrix.increment(timeSeries.trueLabels()[i], timeSeries.predictedLabels()[i]);
            if (timeSeries.trueLabels()[i] == timeSeries.predictedLabels()[i]) {
                truePositives.increment(timeSeries.trueLabels()[i]);
                correctPredictions++;
            } else {
                falseNegatives.increment(timeSeries.trueLabels()[i]);
                falsePositives.increment(timeSeries.predictedLabels()[i]);
            }
        }
    }

    public Map<String, Object> eval() {
        for (Integer label : allLabels) {
            int labelTP = truePositives.count(label);
            int labelFP = falsePositives.count(label);
            int labelFN = falseNegatives.count(label);
            if ((labelTP + labelFP) == 0) {
                neverPredictedLabels.add(label);
            } else {
                precision.put(label, ((double) labelTP) / (labelTP + labelFP));
            }
            if ((labelTP + labelFN) == 0) {
                neverAppearedLabels.add(label);
            } else {
                recall.put(label, ((double) labelTP) / (labelTP + labelFN));
            }
        }
        double totalPrecision = 0;
        double totalRecall = 0;
        for (Double labelPrecision : precision.values()) {
            totalPrecision += labelPrecision;
        }
        for (Double labelRecall : recall.values()) {
            totalRecall += labelRecall;
        }
        final double finalPrecision = totalPrecision / precision.size();
        final double finalRecall = totalRecall / recall.size();
        stats.put("accuracy", (double) correctPredictions / totalPredictions);
        stats.put("precision", finalPrecision);
        stats.put("recall", finalRecall);
        stats.put("f1", (2.0 * finalPrecision * finalRecall) / (finalPrecision + finalRecall));
        stats.put("never_appeared", neverAppearedLabels.toArray());
        stats.put("never_predicted", neverPredictedLabels.toArray());
        stats.putAll(confusionMatrix.toMap());
        return stats;
    }

    public String evalString() {
        assert stats.size() > 0 : "eval() should be called first";
        StringBuilder statsBuilder = new StringBuilder();
        statsBuilder.append(String.format("\t%f", (float) stats.get("accuracy")));
        statsBuilder.append(String.format("\t%f", (float) stats.get("precision")));
        statsBuilder.append(String.format("\t%f", (float) stats.get("recall")));
        statsBuilder.append(String.format("\t%f", (float) stats.get("f1")));
        statsBuilder.append(String.format("\t%s", Arrays.toString((int[]) stats.get("never_appeared"))));
        statsBuilder.append(String.format("\t%s", Arrays.toString((int[]) stats.get("never_predicted"))));
        for (int i : allLabels) {
            for (int j : allLabels) {
                statsBuilder.append(String.format("\t%d",
                        (int) stats.get(String.format("%d_actual_%d_predicted", i, j))));
            }
        }
        return statsBuilder.toString();
    }

    private class Counter {
        private Map<Integer, Integer> backingMap;

        private Counter() {
            backingMap = new HashMap<>();
        }

        private void increment(int label) {
            if (backingMap.containsKey(label)) {
                backingMap.put(label, backingMap.get(label) + 1);
            } else {
                backingMap.put(label, 1);
            }
        }

        private int count(int label) {
            Integer labelCount = backingMap.get(label);
            return labelCount != null ? labelCount : 0;
        }
    }

    private class ConfusionMatrix {
        private Map<Integer, Counter> backingMap;
        private List<Integer> allLabels;

        private ConfusionMatrix(List<Integer> allLabels) {
            backingMap = new HashMap<>();
            this.allLabels = allLabels;
            for (Integer label : allLabels) {
                backingMap.put(label, new Counter());
            }
        }

        private void increment(int trueLabel, int predictedLabel) {
            backingMap.get(trueLabel).increment(predictedLabel);
        }

        private int count(int trueLabel, int predictedLabel) {
            return backingMap.get(trueLabel).count(predictedLabel);
        }

        private Map<String, Integer> toMap() {
            Map<String, Integer> confusionMatrixMap = new HashMap<>();
            for (int i : allLabels) {
                for (int j : allLabels) {
                    confusionMatrixMap.put(String.format("%d_actual_%d_predicted", i, j), count(i, j));
                }
            }
            return confusionMatrixMap;
        }
    }

    public static double estimateFromGraph(ComputationGraph graph, MultiDataSetIterator iterator, int numLabels,
                                           String metricName, int outputIndex, long scoreN) {
        TimeSeriesPerformanceCalculator calculator = new TimeSeriesPerformanceCalculator(numLabels);
        int sequenceCount = 0;
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            if (next.hasMaskArrays()) {
                graph.setLayerMaskArrays(next.getFeaturesMaskArrays(), next.getLabelsMaskArrays());
            }
            INDArray allTrueLabels = graph.output(next.getFeatures())[outputIndex];
            int numExamples = next.getFeatures(outputIndex).size(0);
            for (int sequenceIdx = 0; sequenceIdx < numExamples; sequenceIdx++) {
                INDArray trueLabel = next.getLabels(outputIndex).getRow(sequenceIdx);
                INDArray predictedLabel = allTrueLabels.getRow(sequenceIdx);
                TimeSeriesPrediction prediction = new TimeSeriesPrediction(trueLabel, predictedLabel);
                calculator.addTimeSeries(prediction);
                sequenceCount++;
            }
            if (sequenceCount > scoreN) {
                break;
            }
        }
        return (double) calculator.eval().get(metricName);
    }
}
