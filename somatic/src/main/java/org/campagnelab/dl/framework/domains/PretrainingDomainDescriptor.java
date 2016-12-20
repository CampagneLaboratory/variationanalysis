package org.campagnelab.dl.framework.domains;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.framework.mappers.pretraining.RNNPretrainingFeatureMapper;
import org.campagnelab.dl.framework.mappers.pretraining.RNNPretrainingLabelMapper;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Arrays;
import java.util.function.Function;

/**
 * Domain descriptor specifically in use for pretraining. Modifies a delegate domain
 * descriptor to work with pretraining tasks.
 *
 * Graph specified by delegate domain descriptor should have node numbers (nIn, nOut) set by the
 * getNumInputs, getNumComponents, and getNumOutputs methods, and must only have one
 * input layer
 *
 * Created by joshuacohen on 12/7/16.
 */
public class PretrainingDomainDescriptor<RecordType> extends DomainDescriptor<RecordType> {
    static private final Logger LOG = LoggerFactory.getLogger(PretrainingDomainDescriptor.class);

    private DomainDescriptor<RecordType> delegate;
    private Function<RecordType, Integer> recordToSequenceLength;
    private TrainingArguments args;
    private String inputName;

    /**
     * Create a PretrainingDomainDescriptor from a delegate domain descriptor, recordToSequenceLength
     * function (passed into the mappers created by the PretrainingDomainDescriptor) and training arguments
     * @param delegate Delegate domain descriptor
     * @param recordToSequenceLength Record to sequence length function
     * @param args Training arguments
     */
    public PretrainingDomainDescriptor(DomainDescriptor<RecordType> delegate,
                                       Function<RecordType, Integer> recordToSequenceLength,
                                       TrainingArguments args) {
        this.delegate = delegate;
        this.recordToSequenceLength = recordToSequenceLength;
        this.args = args;
        ComputationGraphAssembler delegateAssembler = delegate.getComputationalGraph();
        this.computationGraphAssembler = delegateAssembler;
        delegate.computationGraphAssembler = delegateAssembler;
        String[] delegateGraphInputs = delegateAssembler.getInputNames();
        if (delegateGraphInputs.length != 1) {
            LOG.warn("Graph should only have one input");
        }
        inputName = delegateGraphInputs[0];
    }

    /**
     * Creates an instance of a RNNPretrainingFeatureMapper from the (required to be 2D) feature mapper
     * from the delegate domain descriptor
     * @param inputName The name of a graph input. Must match an input of the computational graph.
     * @return RNNPretrainingFeatureMapper instance
     */
    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        if (!inputName.equals(this.inputName)) {
            LOG.warn("Invalid input name; given {} but should be {}", inputName, this.inputName);
        }
        FeatureMapper delegateFeatureMapper = delegate.getFeatureMapper(inputName);
        return new RNNPretrainingFeatureMapper<>(delegateFeatureMapper, args.eosIndex, recordToSequenceLength);
    }

    /**
     * Creates an instance of a RNNPretrainingLabelMapper from the (required to be 2D) feature mapper
     * from the delegate domain descriptor (since the labels in pretraining correspond to the features
     * from the regular training task)
     * @param outputName The name of a graph output. Must match an output of the computational graph.
     * @return RNNPretrainingLabelMapper instance
     */
    @Override
    public LabelMapper getLabelMapper(String outputName) {
        FeatureMapper delegateLabelMapper = delegate.getFeatureMapper(inputName);
        return new RNNPretrainingLabelMapper<>(delegateLabelMapper, args.eosIndex, recordToSequenceLength);
    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {
        return delegate.getPredictionInterpreter(outputName);
    }

    @Override
    public Function<String, ? extends Iterable<RecordType>> getRecordIterable() {
        return delegate.getRecordIterable();
    }

    @Override
    public int[] getNumInputs(String inputName) {
        if (!inputName.equals(this.inputName)) {
            LOG.warn("Invalid input name; given {} but should be {}", inputName, this.inputName);
        }
        int[] delegateNumInputs = delegate.getNumInputs(inputName).clone();
        if (delegateNumInputs.length != 2) {
            throw new IllegalArgumentException("Delegate number of inputs should be two dimensional");
        }
        boolean addEosIndex = (args.eosIndex != null && args.eosIndex == delegateNumInputs[0]) || args.eosIndex == null;
        if (addEosIndex) delegateNumInputs[0]++;
        delegateNumInputs[1] *= 2;
        delegateNumInputs[1]++;
        return delegateNumInputs;
    }

    @Override
    public int[] getNumOutputs(String outputName) {
        return getNumInputs(inputName);
    }

    @Override
    public int[] getNumMaskInputs(String inputName) {
        if (!inputName.equals(this.inputName)) {
            LOG.warn("Invalid input name; given {} but should be {}", inputName, this.inputName);
        }
        int[] delegateNumMaskInputs = delegate.getNumMaskInputs(inputName).clone();
        if (delegateNumMaskInputs.length != 1) {
            throw new IllegalArgumentException("Delegate mask should be one dimensional");
        }
        delegateNumMaskInputs[0] *= 2;
        delegateNumMaskInputs[0]++;
        return delegateNumMaskInputs;
    }

    @Override
    public int[] getNumMaskOutputs(String outputName) {
        return getNumMaskInputs(inputName);
    }

    @Override
    public int getNumHiddenNodes(String componentName) {
        return delegate.getNumHiddenNodes(componentName);
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        return delegate.getNumRecords(recordFiles);
    }

    @Override
    public ComputationGraphAssembler getComputationalGraph() {
        return delegate.getComputationalGraph();
    }

    @Override
    public LossFunctions.LossFunction getOutputLoss(String outputName) {
        return LossFunctions.LossFunction.MCXENT;
    }
}
