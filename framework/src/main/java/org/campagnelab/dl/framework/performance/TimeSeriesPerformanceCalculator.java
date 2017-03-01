package org.campagnelab.dl.framework.performance;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.domains.prediction.TimeSeriesPrediction;
import org.campagnelab.dl.framework.domains.prediction.TimeSeriesPredictionInterpreter;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.ImmutablePair;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Calculates stats for time series predictions
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
    private Map<Integer, Double> labelPrecisions;
    private Map<Integer, Double> labelRecalls;
    private double mcPrecision;
    private double mcRecall;
    private double mcAccuracy;
    private double mcF1Score;
    private ConfusionMatrix confusionMatrix;
    private List<Integer> allLabels;
    private boolean evalCalled;


    public TimeSeriesPerformanceCalculator(List<Integer> allLabels) {
        truePositives = new Counter();
        falsePositives = new Counter();
        falseNegatives = new Counter();
        neverPredictedLabels = new HashSet<>();
        neverAppearedLabels = new HashSet<>();
        labelPrecisions = new HashMap<>();
        labelRecalls = new HashMap<>();
        confusionMatrix = new ConfusionMatrix(allLabels);
        totalPredictions = 0;
        this.allLabels = allLabels;
    }

    public TimeSeriesPerformanceCalculator(int numLabels) {
        this(IntStream.range(0, numLabels).boxed().collect(Collectors.toList()));
    }

    public TimeSeriesPerformanceCalculator addTimeSeries(TimeSeriesPrediction timeSeries) {
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
        return this;
    }

    public TimeSeriesPerformanceCalculator eval() {
        for (Integer label : allLabels) {
            int labelTP = truePositives.count(label);
            int labelFP = falsePositives.count(label);
            int labelFN = falseNegatives.count(label);
            if ((labelTP + labelFP) == 0) {
                neverPredictedLabels.add(label);
            } else {
                labelPrecisions.put(label, ((double) labelTP) / (labelTP + labelFP));
            }
            if ((labelTP + labelFN) == 0) {
                neverAppearedLabels.add(label);
            } else {
                labelRecalls.put(label, ((double) labelTP) / (labelTP + labelFN));
            }
        }
        double totalPrecision = 0;
        double totalRecall = 0;
        for (Double labelPrecision : labelPrecisions.values()) {
            totalPrecision += labelPrecision;
        }
        for (Double labelRecall : labelRecalls.values()) {
            totalRecall += labelRecall;
        }
        mcPrecision = totalPrecision / labelPrecisions.size();
        mcRecall = totalRecall / labelRecalls.size();
        mcAccuracy = (double) correctPredictions / totalPredictions;
        mcF1Score = (2.0 * mcPrecision * mcRecall) / (mcPrecision + mcRecall);
        evalCalled = true;
        return this;
    }

    public double getMetric(String metricName) {
        assert evalCalled : "eval() should be called first";
        switch (metricName) {
            case "f1":
                return mcF1Score;
            case "accuracy":
                return mcAccuracy;
            case "precision":
                return mcPrecision;
            case "recall":
                return mcRecall;
            default:
                throw new UnsupportedOperationException("Unknown metric name");
        }
    }

    public Set<Integer> getNeverAppearedOrPredictedSet(String setName) {
        assert evalCalled : "eval() should be called first";
        switch (setName) {
            case "never_appeared":
                return neverAppearedLabels;
            case "never_predicted":
                return neverPredictedLabels;
            default:
                throw new UnsupportedOperationException("Unknown set name");
        }
    }

    public Map<Pair<Integer, Integer>, Integer> getConfusionMatrix() {
        assert evalCalled : "eval() should be called first";
        return confusionMatrix.toMap();
    }

    public int countConfusionMatrix(int trueLabel, int predictedLabel) {
        return confusionMatrix.count(trueLabel, predictedLabel);
    }


    public String evalString() {
        StringBuilder statsBuilder = new StringBuilder();
        statsBuilder.append(String.format("\t%f", getMetric("accuracy")));
        statsBuilder.append(String.format("\t%f", getMetric("precision")));
        statsBuilder.append(String.format("\t%f", getMetric("recall")));
        statsBuilder.append(String.format("\t%f", getMetric("f1")));
        statsBuilder.append(String.format("\t%s", getNeverAppearedOrPredictedSet("never_appeared")));
        statsBuilder.append(String.format("\t%s", getNeverAppearedOrPredictedSet("never_predicted")));
        for (int i : allLabels) {
            for (int j : allLabels) {
                statsBuilder.append(String.format("\t%d", confusionMatrix.count(i, j)));
            }
        }
        return statsBuilder.toString();
    }

    private class Counter {
        private Map<Integer, Integer> backingMap;

        public Counter() {
            backingMap = new HashMap<>();
        }

        public void increment(int label) {
            if (backingMap.containsKey(label)) {
                backingMap.put(label, backingMap.get(label) + 1);
            } else {
                backingMap.put(label, 1);
            }
        }

        public int count(int label) {
            Integer labelCount = backingMap.get(label);
            return labelCount != null ? labelCount : 0;
        }
    }

    private class ConfusionMatrix {
        private Map<Integer, Counter> backingMap;
        private List<Integer> allLabels;

        public ConfusionMatrix(List<Integer> allLabels) {
            backingMap = new HashMap<>();
            this.allLabels = allLabels;
            for (Integer label : allLabels) {
                backingMap.put(label, new Counter());
            }
        }

        public void increment(int trueLabel, int predictedLabel) {
            backingMap.get(trueLabel).increment(predictedLabel);
        }

        public int count(int trueLabel, int predictedLabel) {
            return backingMap.get(trueLabel).count(predictedLabel);
        }

        public Map<Pair<Integer, Integer>, Integer> toMap() {
            Map<Pair<Integer, Integer>, Integer> confusionMatrixMap = new HashMap<>();
            for (int i : allLabels) {
                for (int j : allLabels) {
                    confusionMatrixMap.put(new ImmutablePair<>(i, j), count(i, j));
                }
            }
            return confusionMatrixMap;
        }
    }

    public static double estimateFromGraph(ComputationGraph graph, MultiDataSetIterator iterator, int numLabels,
                                           String metricName, int outputIndex, long scoreN,
                                           DomainDescriptor domainDescriptor) {

        return estimateFromGraphHelper(graph, iterator, numLabels, outputIndex, scoreN, domainDescriptor)
                .getMetric(metricName);
    }

    public static Map<String, Double> estimateFromGraph(ComputationGraph graph, MultiDataSetIterator iterator, int numLabels,
                                             int outputIndex, long scoreN, DomainDescriptor domainDescriptor,
                                             String... metrics) {
        TimeSeriesPerformanceCalculator timeSeriesPerformanceCalculator = estimateFromGraphHelper(graph,
                iterator, numLabels, outputIndex, scoreN, domainDescriptor);
        Map<String, Double> graphMetrics = new HashMap<>();
        for (String metric : metrics) {
            if (!metric.equals("score")) {
                graphMetrics.put(metric, timeSeriesPerformanceCalculator.getMetric(metric));
            }
        }
        return graphMetrics;
    }

    private static TimeSeriesPerformanceCalculator estimateFromGraphHelper(ComputationGraph graph, MultiDataSetIterator iterator, int numLabels,
                                                                           int outputIndex, long scoreN,
                                                                           DomainDescriptor domainDescriptor) {
        PredictionInterpreter predictionInterpreter = domainDescriptor.getPredictionInterpreter(
                domainDescriptor.getComputationalGraph().getOutputNames()[outputIndex]);
        if (!(predictionInterpreter instanceof TimeSeriesPredictionInterpreter)) {
            throw new IllegalArgumentException("TimeSeriesPerformanceCalculator needs a TimeSeriesPredictionInterpreter");
        }
        TimeSeriesPredictionInterpreter timeSeriesPredictionInterpreter = (TimeSeriesPredictionInterpreter) predictionInterpreter;
        TimeSeriesPerformanceCalculator calculator = new TimeSeriesPerformanceCalculator(numLabels);
        int sequenceCount = 0;
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            if (next.hasMaskArrays()) {
                graph.setLayerMaskArrays(next.getFeaturesMaskArrays(), next.getLabelsMaskArrays());
            }
            INDArray allPredictedLabels = graph.output(next.getFeatures())[outputIndex];
            INDArray allTrueLabels = next.getLabels(outputIndex);
            int numExamples = next.getFeatures(outputIndex).size(0);
            for (int sequenceIdx = 0; sequenceIdx < numExamples; sequenceIdx++) {
                TimeSeriesPrediction prediction = timeSeriesPredictionInterpreter.interpret(allTrueLabels,
                        allPredictedLabels, sequenceIdx);
                calculator.addTimeSeries(prediction);
                sequenceCount++;
            }
            if (sequenceCount > scoreN) {
                break;
            }
        }
        iterator.reset();
        return calculator.eval();
    }
}
