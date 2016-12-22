package org.campagnelab.dl.framework.domains;

import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.*;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.nd4j.linalg.lossfunctions.LossFunctions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Properties;
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
    private TrainingArguments args;
    private String inputName;
    private String pretrainingFeatureMapperClassName;
    private String pretrainingLabelMapperClassName;
    private Properties pretrainingProperties;

    /**
     * Create a PretrainingDomainDescriptor from a delegate domain descriptor, recordToSequenceLength
     * function (passed into the mappers created by the PretrainingDomainDescriptor) and training arguments
     * @param delegate Delegate domain descriptor
     * @param args Training arguments
     * @param pretrainingFeatureMapperClassName Default feature class name
     * @param pretrainingLabelMapperClassName Default label mapper class name
     * @param pretrainingProperties pretraining properties
     */
    public PretrainingDomainDescriptor(DomainDescriptor<RecordType> delegate,
                                       TrainingArguments args, String pretrainingFeatureMapperClassName,
                                       String pretrainingLabelMapperClassName,
                                       Properties pretrainingProperties) {
        this.delegate = delegate;
        this.args = args;
        ComputationGraphAssembler delegateAssembler = delegate.getComputationalGraph();
        this.computationGraphAssembler = delegateAssembler;
        delegate.computationGraphAssembler = delegateAssembler;
        String[] delegateGraphInputs = delegateAssembler.getInputNames();
        if (delegateGraphInputs.length != 1) {
            throw new IllegalArgumentException("Graph should only have one input");
        }
        inputName = delegateGraphInputs[0];
        super.loadProperties(delegate.domainProperties, delegate.modelProperties);
        this.pretrainingFeatureMapperClassName = pretrainingFeatureMapperClassName;
        this.pretrainingLabelMapperClassName = pretrainingLabelMapperClassName;
        this.pretrainingProperties = pretrainingProperties;
    }

    /**
     * Creates an instance of the default pretrainingFeatureMapper class and configures it
     * @param inputName name of graph input
     * @return FeatureMapper instance configured with pretrainingProperties
     */
    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        try {
            FeatureMapper featureMapper = (FeatureMapper) Class.forName(pretrainingFeatureMapperClassName)
                    .newInstance();
            ConfigurableFeatureMapper configurableFeatureMapper = (ConfigurableFeatureMapper) featureMapper;
            configurableFeatureMapper.configure(pretrainingProperties);
            return featureMapper;
        } catch (Exception e) {
            throw new IllegalArgumentException("Invalid instance of feature mapper", e);
        }

    }

    /**
     * Creates an instance of the default pretrainingLabelMapper class and configures it
     * @param outputName name of graph output
     * @return LabelMapper instance configured with pretrainingProperties
     */
    @Override
    public LabelMapper getLabelMapper(String outputName) {
        try {
            LabelMapper labelMapper = (LabelMapper) Class.forName(pretrainingLabelMapperClassName)
                    .newInstance();
            ConfigurableLabelMapper configurableLabelMapper = (ConfigurableLabelMapper) labelMapper;
            configurableLabelMapper.configure(pretrainingProperties);
            return labelMapper;
        } catch (Exception e) {
            throw new IllegalArgumentException("Invalid instance of feature mapper", e);
        }
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
