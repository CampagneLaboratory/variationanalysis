package org.campagnelab.dl.framework.performance;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

public abstract class PerformanceMetricDescriptor<RecordType> {
    protected DomainDescriptor<RecordType> domainDescriptor;

    public PerformanceMetricDescriptor(DomainDescriptor<RecordType> domainDescriptor) {
        this.domainDescriptor = domainDescriptor;
    }

    /**
     * Return the names of performance metrics to evaluate for this domain.
     *
     * @return
     */
    public abstract String[] performanceMetrics();

    /**
     * Indicate if a larger value of the metric is better performance (e.g., AUC).
     *
     * @param metricName
     * @return True when a larger numerical value indicates better performance. False otherwise (smaller values are best).
     */
    public abstract boolean largerValueIsBetterPerformance(String metricName);

    /**
     * Estimate performance for a metric and return its value.
     *
     * @param metricName
     * @return
     */
    public abstract double estimateMetric(ComputationGraph graph, String metricName,
                                          MultiDataSetIterator dataSetIterator, long scoreN);

    /**
     * A method that can evaluate several metrics in one pass over a validation set. Metric values are
     * returned in the order of the metrics names provided as arguments.
     * Note that the default implementation is not optimal because it performs several
     * passes over the validation set. A custom implementation can often improve on this.
     *
     * @param graph           model.
     * @param dataSetIterator iterator over validation set.
     * @param scoreN          number of examples to score.
     * @param metrics         name of metrics to estimate.
     * @return an array of metric values.
     */
    public double[] estimateMetric(ComputationGraph graph,
                                   MultiDataSetIterator dataSetIterator, long scoreN, String... metrics) {
        DoubleList results = new DoubleArrayList();
        for (String metricName : metrics) {
            results.add(estimateMetric(graph, metricName, dataSetIterator, scoreN));
        }
        return results.toDoubleArray();
    }

    /**
     * Return the name of the metric to use for early stopping.
     *
     * @return
     */
    public abstract String earlyStoppingMetric();
}