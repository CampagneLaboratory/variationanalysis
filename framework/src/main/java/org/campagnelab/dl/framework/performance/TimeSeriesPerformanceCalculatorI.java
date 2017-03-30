package org.campagnelab.dl.framework.performance;

import org.campagnelab.dl.framework.domains.prediction.TimeSeriesPrediction;

/**
 * Created by fac2003 on 3/29/17.
 */
public interface TimeSeriesPerformanceCalculatorI {
    /** Evaluate statistics. Call this method before getMetric(). */
    public TimeSeriesPerformanceCalculatorI eval();

    /**
     * Obtain a performance metric value. This method can be called after eval() has been invoked.
     * @param metricName name of the metric.
     * @return the value of the performance metric.
     */
    public double getMetric(String metricName);

    /**
     * Observe a time series. Must be called at least once before calling eval.
     * @param timeSeries timeSeries prediction.
     * @return
     */
    public void addTimeSeries(TimeSeriesPrediction timeSeries);
}
