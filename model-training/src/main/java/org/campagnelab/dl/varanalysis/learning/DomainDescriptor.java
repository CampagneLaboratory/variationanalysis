package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.glassfish.jersey.internal.util.Producer;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.util.function.Function;

/**
 * A domain descriptor provides information about a modelling domain. Implementations of the domain descriptor
 * define precisely the variable parts of a domain (such as feature and label mappers, loss functions for each
 * label, etc.).
 * Implementing this class provides the information needed to start training models for a new domain. See TrainModel
 * and the TrainModelS example where we define a somatic mutation modeling domain.
 *
 * @param <RecordType> The type of record (i.e., record in .sbi file for instance) that is processed in the domain.
 */
public abstract class DomainDescriptor<RecordType> {
    /**
     * Get the feature mapper for a given model graph input.
     * @param inputName The name of a graph input. Must match an input of the computational graph.
     * @return A feature mapper.
     */
    public abstract FeatureMapper getFeatureMapper(String inputName);
    /**
     * Get the label mapper for a given model graph output.
     * @param outputName The name of a graph output. Must match an output of the computational graph.
     * @return A label mapper.
     */
    public abstract LabelMapper getLabelMapper(String outputName);

    /**
     * Returns a function that converts an input filename to an iterable over records in the file.
     *
     * @return
     */
    public abstract Function<String, ? extends Iterable<RecordType>> getRecordIterable();

    /**
     * Return a computational graph assembler. The assembler can build a computational graph ready for training.
     *
     * @return ComputationalGraphAssembler
     */
    public abstract ComputationalGraphAssembler getComputationalGraph();

    /**
     * Return the dimensions of an input to the graph. If the graph has one input with 10 features, this method should return new int[]{10}.
     * If the graph has one input with two dimensions width=110 and height=120, this method should return new int[]{110,120}.
     * @param inputName input of the computational graph.
     * @return the number of dimensions for each record. 1d records will have one dimension, 2-d records will have 2, etc.
     */
    public abstract int[] getNumInputs(String inputName);
    /**
     * Return the dimensions of an output of the graph.
     * @param outputName output of the computational graph.
     * @return the number of dimensions for each record. 1d records will have one dimension, 2-d records will have 2, etc.
     */
    public abstract int[] getNumOutputs(String outputName);
    /**
     * Return the number of hidden nodes for a component of the graph. The number will be used to configure the graph.
     * The number may vary
     * @param componentName component of the computational graph.
     * @return the number of hidden nodes to allocate to the component.
     */
    public abstract int getNumHiddenNodes(String componentName);

    /**
     * Return the output loss for an output.
     * @param outputName
     * @return
     */
    public abstract LossFunctions.LossFunction getOutputLoss(String outputName);

    /**
     * Return the number of records across the record files (i.e., .sbi files).
     *
     * @return number of records.
     */
    public abstract long getNumRecords(String[] recordFiles);

    // The following provide default implementations suitable when training with only Loss score.
    public PerformanceMetricDescriptor<RecordType> performanceDescritor() {
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
                return estimateScore(graph, metricName, dataSetIterator, scoreN);
            }


            @Override
            public String earlyStoppingMetric() {
                return "score";
            }
        };
    }

    protected double estimateScore(ComputationGraph graph, String metricName, MultiDataSetIterator dataSetIterator, long scoreN) {
        switch (metricName) {
            case "score":

                double score = 0;
                long nBatch = 0;
                long nExamples = 0;
                while (dataSetIterator.hasNext()) {
                    MultiDataSet ds = dataSetIterator.next();
                    score += graph.score(ds);
                    nBatch += 1;
                    nExamples += ds.getFeatures()[0].size(0);
                    if (nExamples > scoreN) break;
                }
                return score / nBatch;
            default:
                throw new IllegalArgumentException("metric name not recognized: " + metricName);
        }
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