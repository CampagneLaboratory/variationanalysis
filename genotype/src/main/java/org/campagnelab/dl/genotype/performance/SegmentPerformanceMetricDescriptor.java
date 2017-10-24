package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

/**
 * Estimate performance of genotype calling over segments. Used when training with .ssi/.ssip files.
 */
public class SegmentPerformanceMetricDescriptor extends PerformanceMetricDescriptor<SegmentInformationRecords.SegmentInformation> {
    public SegmentPerformanceMetricDescriptor(DomainDescriptor<SegmentInformationRecords.SegmentInformation> domainDescriptor) {
        super(domainDescriptor);
        helper = new SegmentTrainingPerformanceHelper(domainDescriptor);
    }

    SegmentTrainingPerformanceHelper helper;

    @Override
    public String[] performanceMetrics() {

        ObjectArrayList<String> metrics = new ObjectArrayList<>();
        metrics.add("validationScore");
        metrics.addAll(helper.getMetricNames());
        return metrics.toArray(new String[metrics.size()]);
    }

    @Override
    public boolean largerValueIsBetterPerformance(String metricName) {

        switch (metricName) {
            case "validationScore":
                return false;
            case "accuracy":
                return true;
            default:
                System.out.println("ordering not set for metric: "+metricName);
                return false;

        }
    }

    @Override
    public double estimateMetric(ComputationGraph graph, String metricName, MultiDataSetIterator dataSetIterator, long scoreN) {

        return 0;
    }

    @Override
    public String earlyStoppingMetric() {
        return "validationScore";
    }

    public double[] estimateMetric(ComputationGraph graph,
                                   MultiDataSetIterator dataSetIterator, long scoreN, String... metrics) {

        helper.estimateWithGraph(dataSetIterator, graph,
                index -> index > scoreN
                            /* first output represents probabilityIsCalled of mutation */);
        return helper.getMetricValues();
    }

}
