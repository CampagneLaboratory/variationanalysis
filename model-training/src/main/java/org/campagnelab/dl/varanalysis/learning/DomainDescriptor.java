package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.campagnelab.dl.varanalysis.learning.iterators.MultiDataSetRecordIterator;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.glassfish.jersey.internal.util.Producer;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.util.function.Function;

public abstract class DomainDescriptor<RecordType> {

    public abstract FeatureMapper getFeatureMapper(String inputName);

    public abstract LabelMapper getLabelMapper(String outputName);

    /**
     * Returns a function that converts an input filename to an iterable over records in the file.
     *
     * @return
     */
    public abstract Function<String, ? extends Iterable<RecordType>> getRecordIterable();

    /**
     * Return a computational graph assembler.
     *
     * @return
     */
    public abstract ComputationalGraphAssembler getComputationalGraph();

    public abstract int[] getNumInputs(String inputName);

    public abstract int[] getNumOutputs(String outputName);

    public abstract int getNumHiddenNodes(String componentName);

    public abstract LossFunctions.LossFunction getOutputLoss(String outputName);

    public  PerformanceMetricDescriptor<RecordType> performanceDescritor() {
        return new PerformanceMetricDescriptor<RecordType>() {
            @Override
            public String[] performanceMetrics() {
                return new String[]{"score"};

            }

            @Override
            public boolean largerValueIsBetterPerformance(String metricName) {
                return false;
            }

            @Override
            public double estimateMetric(ComputationGraph graph, String metricName, MultiDataSetIterator dataSetIterator, long scoreN) {
                double score=0;
                long nBatch=0;
                long nExamples=0;
                while (dataSetIterator.hasNext()) {
                    MultiDataSet ds = dataSetIterator.next();
                    score+=graph.score(ds);
                    nBatch+=1;
                    nExamples+=ds.getFeatures()[0].size(0);
                    if (nExamples>scoreN) break;
                }
                return score/nBatch;

            }


            @Override
            public String earlyStoppingMetric() {
                return "score";
            }
        };
    }

    public int[] getInputShape(int size, String inputName) {
        return getShape(size, () -> getNumInputs(inputName));
    }

    public int[] getLabelShape(int size, String outputName) {
        return getShape(size, () -> getNumOutputs(outputName));
    }

    public int[] getShape(int size, Producer<int[]> p) {
        int[] numInputs = p.call();
        assert numInputs.length <= 2;
        switch (numInputs.length) {
            case 1:
                return new int[]{size, numInputs[0]};
            case 2:
                return new int[]{size, numInputs[0], numInputs[1]};
            default:
                throw new UnsupportedOperationException();
        }
    }


}