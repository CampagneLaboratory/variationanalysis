package org.campagnelab.dl.framework.performance;

import it.unimi.dsi.fastutil.objects.Object2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

/**
 * Write performance logs under the model directory. Supports multiple models produced during training,
 * where models are identified by a prefix: best, bestAUC, final, 1-.
 * Created by fac2003 on 7/24/16.
 */
public class PerformanceLogger {
    private static final String perfFilenameFormat = "%s-perf-log.tsv";
    private String directory;
    private String conditionId;
    private double[] bestPerformances;
    private boolean[] performanceLargeIsBest;
    private String[] performanceNames;

    /**
     * Define performance metrics where the first element of a tuple is the metric name, the second
     * a boolean that is true when large values of the metric are better than small ones.
     *
     * @param metrics A set of tuples describing metrics to collect.
     */
    public void definePerformances(Metric... metrics) {
        performanceLargeIsBest = new boolean[metrics.length];
        performanceNames = new String[metrics.length];
        bestPerformances = new double[metrics.length];
        int index = 0;
        for (Metric metric : metrics) {
            performanceNames[index] = metric.name;
            performanceLargeIsBest[index] = metric.largerIsBetter;
            bestPerformances[index] = performanceLargeIsBest[index] ? Double.NEGATIVE_INFINITY : Double.MAX_VALUE;
            index += 1;
        }
    }

    private double bestScore = Float.MAX_VALUE;
    private double bestAUC = -1;

    /**
     * Return the best score seen so far. Smaller score is best.
     *
     * @return
     */
    public double getBestScore() {
        int index = 0;
        for (String metricName : performanceNames) {
            if ("score".equals(metricName)) {
                return bestPerformances[index];
            }
            index++;
        }
        return Double.NaN;
    }

    /**
     * Return the best AUC seen so far. Larger  AUC is best.
     *
     * @return
     */
    public double getBestAUC() {
        return bestAUC;
    }

    /**
     * Returns the best value of the performance metric.
     *
     * @param performanceName Name of the performance metric.
     * @return best performance seen.
     */
    public double getBest(String performanceName) {
        int index = 0;
        for (String name : performanceNames) {
            if (name.equals(performanceName)) {
                return bestPerformances[index];
            }
            index += 1;
        }
        throw new IllegalArgumentException("The metric name was not defined: " + performanceName);
    }

    /**
     * @param directory Directory to save model performance logs to.
     */
    public PerformanceLogger(String directory) {
        this.directory = directory;

        final File dir = new File(directory);
        if (!dir.exists() && !dir.mkdirs()) {
            throw new IllegalArgumentException("Unable to create log directory at " + directory + " current directory is " + new File(".").getAbsolutePath());
        }
    }

    public void setCondition(String conditionId) {
        this.conditionId = conditionId;
    }

    public void clear() {
        log.clear();
    }

    /**
     * Log the performance of a model. The model is identified by a prefix. Use "best" if you are only
     * tracking the performance of the best model, by score. Use bestAUC if you are tracking the performance
     * of the model that got best AUC on the test set. Use final for the model created at the end of the
     * training process, irrespective of performance.
     *
     * @param prefix          Identifies a specific model.
     * @param numExamplesUsed The number of training examples seen by the model. Note that examples used in training several times count several times.
     * @param epoch           The number of epochs used to train the model.
     * @param metricValues    values of performance metrics.
     */
    public void logMetrics(String prefix, long numExamplesUsed, int epoch, double... metricValues) {
        ObjectArrayList<Performance> defaultValue = new ObjectArrayList<>();
        for (int metricIndex = 0; metricIndex < metricValues.length; metricIndex++) {
            if (metricValues[metricIndex] == metricValues[metricIndex]) {
                // not NaN
                if (performanceLargeIsBest[metricIndex]) {
                    bestPerformances[metricIndex] = Math.max(bestPerformances[metricIndex], metricValues[metricIndex]);
                } else {
                    bestPerformances[metricIndex] = Math.min(bestPerformances[metricIndex], metricValues[metricIndex]);
                }
            }
        }
        log.getOrDefault(prefix, defaultValue).add(new Performance(numExamplesUsed, epoch, performanceNames, metricValues));
        if (defaultValue.size() > 0) {
            log.put(prefix, defaultValue);
        }
    }

    public void logTrainingScore(String prefix, int epoch, double trainingScore) {
        List<Performance> values = log.get(prefix);
        Performance perf = values.stream().filter(performance -> performance.epoch == epoch).findFirst().get();
        perf.trainingScore = trainingScore;

    }

    public void show(String prefix) {
        List<Performance> perfs = log.get(prefix);
        Performance perf = perfs.get(perfs.size() - 1);
        System.out.printf("%d\t%f\t%s%n", perf.epoch, perf.trainingScore, perf.formatValues());
    }

    /**
     * Log the performance of a model. The model is identified by a prefix. Use "best" if you are only
     * tracking the performance of the best model, by score. Use bestAUC if you are tracking the performance
     * of the model that got best AUC on the test set. Use final for the model created at the end of the
     * training process, irrespective of performance.
     *
     * @param prefix          Identifies a specific model.
     * @param numExamplesUsed The number of training examples seen by the model. Note that examples used in training several times count several times.
     * @param epoch           The number of epochs used to train the model.
     * @param score           The  score obtained at numExamplesUsed and epoch for the model.
     * @param auc             The AUC on the test set, or NaN.
     */
    public void log(String prefix, long numExamplesUsed, int epoch, double score, double auc) {
        ObjectArrayList<Performance> defaultValue = new ObjectArrayList<>();
        bestScore = Math.min(bestScore, score);
        bestAUC = Math.max(bestAUC, auc);
        log.getOrDefault(prefix, defaultValue).add(new Performance(numExamplesUsed, epoch, score, auc));
        if (defaultValue.size() > 0) {
            log.put(prefix, defaultValue);
        }
    }

    /**
     * Return the epoch when the best performance was obtained. The maximum epoch over all recorded performance logs
     * is returned. (Do not log performance when the performance is not best so far if you use this method.)
     */
    public int getBestEpoch(String prefix) {
        int epoch = -1;
        List<Performance> performances = log.get(prefix);
        if (performances == null) {
            return -1;
        }
        for (Performance per : performances) {
            epoch = Math.max(epoch, per.epoch);
        }
        return epoch;
    }

    /**
     * Log the performance of a model. The model is identified by a prefix. Use "best" if you are only
     * tracking the performance of the best model, by score. Use bestAUC if you are tracking the performance
     * of the model that got best AUC on the test set. Use final for the model created at the end of the
     * training process, irrespective of performance.
     *
     * @param prefix          Identifies a specific model.
     * @param numExamplesUsed The number of training examples used to train the model so far. Note that
     *                        reused examples are counted again.
     * @param epoch           The number of epochs used to train the model.
     * @param score           The  score obtained at numExamplesUsed and epoch for the model.
     */
    public void log(String prefix, long numExamplesUsed, int epoch, double score) {
        log.getOrDefault(prefix, new ObjectArrayList<>()).add(new Performance(numExamplesUsed, epoch, score, -1));
    }

    /**
     * Write performance log for all model prefixes.
     *
     * @throws IOException
     */
    public void write() throws IOException {
        for (String prefix : log.keySet()) {
            write(prefix);
        }
    }

    /**
     * Write performance log for one model prefix.
     *
     * @throws IOException
     */
    public void write(String prefix) throws IOException {

        Writer writer = new FileWriter(directory + "/" +
                String.format(PerformanceLogger.perfFilenameFormat, prefix));
        try {
            writeHeaders(writer);
            List<Performance> perfs = log.get(prefix);
            if (perfs == null) {
                return;
            }
            for (Performance perf : perfs) {

                writer.write(String.format("%d\t%d\t%f\t%s",
                        perf.numExamplesUsed, perf.epoch, perf.trainingScore, perf.formatValues()));

                if (conditionId != null) {
                    writer.write("\t" + conditionId);
                }
                writer.write("\n");
            }
            writer.flush();
        } finally {
            writer.close();
        }

    }

    private void writeHeaders(Writer writer) throws IOException {
        writer.write("numExamplesUsed\tepoch\ttrainingScore\t" + getMetricHeader());
        if (conditionId != null) {
            writer.write("\tcondition");
        }
        writer.write("\n");
    }


    private Object2ObjectMap<String, List<Performance>> log = new Object2ObjectAVLTreeMap<>();

    public String getMetricHeader() {
        if (performanceNames == null) {
            return "score\tAUC";
        } else {
            String result = "epoch\ttrainingScore\t";
            int index = 0;
            for (String name : performanceNames) {
                result += name;
                index += 1;
                if (index < performanceNames.length) {
                    result += "\t";
                }
            }
            return result;
        }
    }

    private class Performance {

        long numExamplesUsed;
        int epoch;
        double[] performanceValues;
        String[] performanceNames;
        /**
         * The score obtained on the training set.
         */
        public double trainingScore;

        public Performance(long numExamplesUsed, int epoch, double score, double auc) {
            this.numExamplesUsed = numExamplesUsed;
            this.epoch = epoch;
            this.score = score;
            this.auc = auc;
        }

        public Performance(long numExamplesUsed, int epoch, String[] performanceNames, double... performanceValues) {
            this.numExamplesUsed = numExamplesUsed;
            this.epoch = epoch;
            this.performanceValues = performanceValues.clone();
            this.performanceNames = performanceNames;
        }

        double score;
        double auc;

        public String formatValues() {
            if (performanceNames == null) {
                return String.format("%f\t%f", score, auc);
            } else {
                String result = "";
                int index = 0;
                for (double value : performanceValues) {
                    if (Math.round(value) == value) {
                        // it's an integer, print it as such:
                        result += String.format("%d", (int)Math.round(value));
                    } else {

                        result += String.format("%f", value);
                    }
                    index += 1;
                    if (index < performanceValues.length) {
                        result += "\t";
                    }
                }
                return result;
            }
        }
    }
}
