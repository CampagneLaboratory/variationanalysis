package org.campagnelab.dl.varanalysis.learning.models;

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

    /**
     * @param directory Directory to save model performance logs to.
     */
    public PerformanceLogger(String directory) {
        this.directory = directory;
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
     * @param prefix    Identifies a specific model.
     * @param iteration The iteration seen by the model.
     * @param epoch     The number of epochs used to train the model.
     * @param score     The  score obtained at iteration and epoch for the model.
     * @param auc       The AUC on the test set, or NaN.
     */
    public void log(String prefix, int iteration, int epoch, double score, double auc) {
        ObjectArrayList<Performance> defaultValue = new ObjectArrayList<>();
        log.getOrDefault(prefix, defaultValue).add(new Performance(iteration, epoch, score, auc));
        if (defaultValue.size() > 0) {
            log.put(prefix, defaultValue);
        }
    }

    /**
     * Log the performance of a model. The model is identified by a prefix. Use "best" if you are only
     * tracking the performance of the best model, by score. Use bestAUC if you are tracking the performance
     * of the model that got best AUC on the test set. Use final for the model created at the end of the
     * training process, irrespective of performance.
     *
     * @param prefix    Identifies a specific model.
     * @param iteration The iteration seen by the model.
     * @param epoch     The number of epochs used to train the model.
     * @param score     The  score obtained at iteration and epoch for the model.
     */
    public void log(String prefix, int iteration, int epoch, double score) {
        log.getOrDefault(prefix, new ObjectArrayList<>()).add(new Performance(iteration, epoch, score, Float.NaN));
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
        new File(directory).mkdirs();
        Writer writer = new FileWriter(directory + "/" +
                String.format(PerformanceLogger.perfFilenameFormat, prefix));
        try {
            writeHeaders(writer);
            List<Performance> perfs = log.get(prefix);
            if (perfs == null) {
                return;
            }
            for (Performance perf : perfs) {
                writer.write(String.format("%d\t%d\t%f\t%f",
                        perf.iteration, perf.epoch, perf.score, perf.auc));
                if (conditionId != null) {
                    writer.write("\t"+conditionId);
                }
                writer.write("\n");
            }
            writer.flush();
        } finally {
            writer.close();
        }

    }

    private void writeHeaders(Writer writer) throws IOException {
        writer.write("iteration\tepoch\tscore\tAUC");
        if (conditionId!=null) {
            writer.write("\tcondition");
        }
        writer.write("\n");
    }


    private Object2ObjectMap<String, List<Performance>> log = new Object2ObjectAVLTreeMap<>();

    private class Performance {

        int iteration;
        int epoch;

        public Performance(int iteration, int epoch, double score, double auc) {
            this.iteration = iteration;
            this.epoch = epoch;
            this.score = score;
            this.auc = auc;
        }

        double score;
        double auc;
    }
}
