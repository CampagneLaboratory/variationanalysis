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
public abstract class PretrainingDomainDescriptor<RecordType> extends DomainDescriptor<RecordType> {
    static private final Logger LOG = LoggerFactory.getLogger(PretrainingDomainDescriptor.class);

    private DomainDescriptor<RecordType> delegate;
    private String inputName;
    private String pretrainingFeatureMapperClassName;
    private String pretrainingLabelMapperClassName;
    private Integer eosIndex;

    /**
     * Create a PretrainingDomainDescriptor from a pretrained model.
     * Initializes an instance of the delegate.domain_descriptor config property using the same modelPath
     * argument; thus, delegate.domain_descriptor must have a one-argument constructor with a String
     * modelPath as its argument
     * @param modelPath path to pretrained model
     */
    public PretrainingDomainDescriptor(String modelPath) {
        super.loadProperties(modelPath);
        String delegateClassName = modelProperties.getProperty("delegate.domain_descriptor");
        if (delegateClassName == null) {
            throw new RuntimeException("Delegate domain descriptor not set. Should be set by PretrainModel run");
        }
        try {
            Object delegateObj = Class.forName(delegateClassName).getConstructor(String.class).newInstance(modelPath);
            delegate = (DomainDescriptor<RecordType>) delegateObj;
        } catch (Exception e) {
            throw new RuntimeException("Invalid instance of delegate domain descriptor. Should have a constructor" +
                    "from a model path and be parameterized on the same RecordType as the " +
                    "PretrainingDomainDescriptor", e);
        }
        initialize();
    }

    /**
     * Create a PretrainingDomainDescriptor from a delegate domain descriptor
     * @param delegate Delegate domain descriptor
     */
    public PretrainingDomainDescriptor(DomainDescriptor<RecordType> delegate) {
        this.delegate = delegate;
        initialize();
    }

    /**
     * Initializes fields from peroperties, abstract method calls
     * Method called by both superclass constructors
     */
    private void initialize() {
        Properties pretrainModelProps = pretrainingModelProperties();
        if (pretrainModelProps == null) pretrainModelProps = new Properties();
        Properties pretrainDomainProps = pretrainingDomainProperties();
        if (pretrainDomainProps == null) pretrainDomainProps = new Properties();
        String eosStringFromModel = delegate.modelProperties.getProperty("delegate.eos_index");
        if (eosStringFromModel != null) modelProperties.setProperty("delegate.eos_index", eosStringFromModel);
        String architectureFromDomain = delegate.domainProperties.getProperty("net.architecture.classname");
        if (architectureFromDomain != null) domainProperties.setProperty("net.architecture.classname",
                architectureFromDomain);
        modelProperties.putAll(pretrainModelProps);
        domainProperties.putAll(pretrainDomainProps);
        String eosString = modelProperties.getProperty("delegate.eos_index");
        eosIndex = eosString.equals("null") ? null : Integer.parseInt(eosString);
        pretrainingFeatureMapperClassName = pretrainingFeatureMapperClassName();
        pretrainingLabelMapperClassName = pretrainingLabelMapperClassName();
        ComputationGraphAssembler delegateAssembler = delegate.getComputationalGraph();
        String[] delegateGraphInputs = delegateAssembler.getInputNames();
        if (delegateGraphInputs.length != 1) {
            throw new IllegalArgumentException("Graph should only have one input");
        }
        inputName = delegateGraphInputs[0];
        this.computationGraphAssembler = delegateAssembler;
        delegate.computationGraphAssembler = delegateAssembler;
        initializeArchitecture();
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
            configurableFeatureMapper.configure(modelProperties);
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
            configurableLabelMapper.configure(modelProperties);
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
        boolean addEosIndex = (eosIndex != null && eosIndex == delegateNumInputs[0]) || eosIndex == null;
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

    /**
     * Fully qualified name of default feature mapper for pretraining. Must have a zero-argument constructor
     * and implement ConfigurableFeatureMapper
     * @return default pretraining feature mapper
     */
    public abstract String pretrainingFeatureMapperClassName();

    /**
     * Fully qualified name of default label mapper for pretraining. Must have a zero-argument constructor
     * and implement ConfigurableLabelMapper
     * @return default pretraining label mapper
     */
    public abstract String pretrainingLabelMapperClassName();

    /**
     * Properties object with any relevant model properties needed to create pretraining feature and label mappers
     * Should contain delegate.eos_index if not defined already in domain
     * @return properties object
     */
    public abstract Properties pretrainingModelProperties();

    /**
     * Properties object with any relevant properties needed to create pretraining feature and label mappers
     * Should contain net.architecture.classname if not defined already in domain
     * @return properties object
     */
    public abstract Properties pretrainingDomainProperties();
}
