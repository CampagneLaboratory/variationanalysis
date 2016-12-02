package org.campagnelab.dl.framework.domains;

import com.google.common.collect.Iterables;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.architecture.graphs.ComputationalGraphAssembler;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Properties;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;

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
     *
     * @param inputName The name of a graph input. Must match an input of the computational graph.
     * @return A feature mapper.
     */
    public abstract FeatureMapper getFeatureMapper(String inputName);

    /**
     * Get the label mapper for a given model graph output.
     *
     * @param outputName The name of a graph output. Must match an output of the computational graph.
     * @return A label mapper.
     */
    public abstract LabelMapper getLabelMapper(String outputName);

    /**
     * Get the prediction/model output interpreter. A prediction interpreter converts the raw
     * INDArray numeric predictions to instances of the BinaryClassPrediction class, in the process converting
     * predictions to a human interpretable representation.
     *
     * @param outputName The name of a graph output. Must match an output of the computational graph.
     * @return A prediction interpreter suitable to interpret the model output predictions.
     */
    public abstract PredictionInterpreter getPredictionInterpreter(String outputName);

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
     *
     * @param inputName input of the computational graph.
     * @return the number of dimensions for each record. 1d records will have one dimension, 2-d records will have 2, etc.
     */
    public abstract int[] getNumInputs(String inputName);

    /**
     * Return the dimensions of an output of the graph.
     *
     * @param outputName output of the computational graph.
     * @return the number of dimensions for each record. 1d records will have one dimension, 2-d records will have 2, etc.
     */
    public abstract int[] getNumOutputs(String outputName);

    /**
     * Return the number of inputs to a computational graph for masking.
     *
     * @param inputName input of the computational graph
     * @return the number of input elements that need masking. Should be a 1D array.
     */
    public abstract int[] getNumMaskInputs(String inputName);

    /**
     * Return the number of outputs from a computational graph for masking.
     *
     * @param outputName output of the computational graph
     * @return the number of output elements that need masking. Should be a 1D array.
     */
    public abstract int[] getNumMaskOutputs(String outputName);

    /**
     * Return the number of hidden nodes for a component of the graph. The number will be used to configure the graph.
     * The number may vary
     *
     * @param componentName component of the computational graph.
     * @return the number of hidden nodes to allocate to the component.
     */
    public abstract int getNumHiddenNodes(String componentName);

    /**
     * Return the output loss for an output.
     *
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
                    double dsSScore = graph.score(ds);
                    if (dsSScore == dsSScore) {
                        // not a NaN
                        score += dsSScore;
                    }
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

    public int[] getInputMaskShape(int size, String inputName) {
        return getShape(size, () -> getNumMaskInputs(inputName));
    }

    public int[] getLabelMaskShape(int size, String outputName) {
        return getShape(size, () -> getNumMaskOutputs(outputName));
    }

    public int[] getShape(int size, Supplier<int[]> p) {
        int[] numInputs = p.get();
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


    public int getNumModelOutputs() {
        return getComputationalGraph().getOutputNames().length;
    }

    public boolean hasOutput(String outputName) {
        for (String name : getComputationalGraph().getOutputNames()) {
            if (outputName.equals(name)) {
                return true;
            }
        }
        return false;
    }

    public FeatureMapper[] featureMappers() {
        FeatureMapper[] mappers = new FeatureMapper[getNumModelInputs()];
        int i = 0;
        for (String inputName : getComputationalGraph().getInputNames()) {
            mappers[i++] = getFeatureMapper(inputName);
        }
        return mappers;
    }

    public int getNumModelInputs() {
        return getComputationalGraph().getInputNames().length;
    }

    protected Properties domainProperties;
    protected Properties modelProperties;

    public void loadProperties(String modelPath) {
        domainProperties = new Properties();
        modelProperties = new Properties();
        String domainPropFilename = ModelLoader.getModelPath(modelPath) + "/domain.properties";
        String modelPropFilename = ModelLoader.getModelPath(modelPath) + "/config.properties";
        try {
            domainProperties.load(new FileReader(domainPropFilename));
            modelProperties.load(new FileReader(modelPropFilename));
        } catch (IOException e) {
            throw new RuntimeException("Unable to load domain properties in model path " + modelPath, e);
        }
    }

    public void loadProperties(Properties domainProperties, Properties sbiProperties) {
        this.domainProperties = new Properties();
        this.modelProperties = new Properties();
        this.domainProperties.putAll(domainProperties);
        this.modelProperties.putAll(sbiProperties);
    }

    public void writeProperties(String modelPath) {
        Properties props = new Properties();
        String propFilename = ModelLoader.getModelPath(modelPath) + "/domain.properties";
        putProperties(props);
        try {
            props.store(new FileWriter(propFilename), "Domain properties created with " + this.getClass().getCanonicalName());
        } catch (IOException e) {
            throw new RuntimeException("Unable to write domain descriptor properties to " + propFilename, e);
        }
    }

    public void putProperties(Properties props) {

        String inputNames[] = getComputationalGraph().getInputNames();
        String outputNames[] = getComputationalGraph().getOutputNames();
        for (String inputName : inputNames) {
            props.put(inputName + ".featureMapper", getFeatureMapper(inputName).getClass().getCanonicalName());
        }
        for (String outputName : outputNames) {
            props.put(outputName + ".labelMapper", getLabelMapper(outputName).getClass().getCanonicalName());
            props.put(outputName + ".predictionInterpreter", getPredictionInterpreter(outputName).getClass().getCanonicalName());
        }

    }

    public Iterable<RecordType> getRecordIterable(List<String> sbiFilenames, int maxRecords) {
        Iterable<RecordType> inputIterable = Iterables.concat(
                sbiFilenames.stream().map(
                        filename -> getRecordIterable().apply(filename)).collect(
                        Collectors.toList()));
        return Iterables.limit(inputIterable, maxRecords);
    }

}