package org.campagnelab.dl.varanalysis.learning.domains;

import org.campagnelab.dl.varanalysis.learning.iterators.MultiDataSetRecordIterator;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

public abstract class PerformanceMetricDescriptor<RecordType>{
    /**
     * Return the names of performance metrics to evaluate for this domain.
     * @return
     */
    public abstract String[]performanceMetrics();

    /**
     * Indicate if a larger value of the metric is better performance (e.g., AUC).
     * @param metricName
     * @return True when a larger numerical value indicates better performance. False otherwise (smaller values are best).
     */
    public abstract boolean largerValueIsBetterPerformance(String metricName);
    /**
     * Estimate performance for a metric and return its value.
     * @param metricName
     * @return
     */
    public abstract double estimateMetric(ComputationGraph graph, String metricName,
                                          MultiDataSetIterator dataSetIterator, long scoreN) ;

        /**
         * Return the name of the metric to use for early stopping.
         * @return
         */
    public abstract String earlyStoppingMetric();
}