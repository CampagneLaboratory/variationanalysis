package org.campagnelab.dl.framework.performance;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by joshuacohen on 11/28/16.
 */
public class TimeSeriesPerformanceCalculator<PredictionType> {
    private Counter truePositives;
    private Counter falsePositives;
    private Counter falseNegatives;
    private Set<PredictionType> encounteredLabels;


    public TimeSeriesPerformanceCalculator() {
        truePositives = new Counter();
        falsePositives = new Counter();
        falseNegatives = new Counter();
        encounteredLabels = new HashSet<>();
    }

    public void addTimeSeries(PredictionType[] trueLabels, PredictionType[] predictedLabels) {
        assert trueLabels.length == predictedLabels.length : "True and predicted labels must have same length";
        for (int i = 0; i < trueLabels.length; i++) {
            encounteredLabels.add(trueLabels[i]);
            encounteredLabels.add(predictedLabels[i]);
            if (trueLabels[i].equals(predictedLabels[i])) {
                truePositives.increment(trueLabels[i]);
            } else {
                falseNegatives.increment(trueLabels[i]);
                falsePositives.increment(predictedLabels[i]);
            }
        }
    }

    private class Counter {
        private Map<PredictionType, Integer> backingMap;
        private Counter() {
            backingMap = new HashMap<>();
        }
        private void increment(PredictionType key) {
            if (backingMap.containsKey(key)) {
                backingMap.put(key, backingMap.get(key) + 1);
            } else {
                backingMap.put(key, 1);
            }
        }
        private Integer get(PredictionType key) {
            return backingMap.get(key);
        }
    }
}
