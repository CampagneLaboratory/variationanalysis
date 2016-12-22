package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.PretrainingDomainDescriptor;

import java.io.IOException;
import java.util.Properties;
import java.util.function.Function;

/**
 * A tool that can be used to create models to perform pretraining on. Creates the domain
 * description using a regular TrainModel passed in as a delegate.
 *
 * Requires implementation of the recordToSequenceLength and createPretrainingDomainDescriptor functions
 *
 * Created by joshuacohen on 12/7/16.
 */
public abstract class PretrainModel<RecordType> extends TrainModel<RecordType> {
    private TrainModel<RecordType> delegate;

    /**
     * Create a pretraining model by using another TrainModel as a delegate
     * @param delegate delegate TrainModel
     */
    public PretrainModel(TrainModel<RecordType> delegate) {
        this.delegate = delegate;
    }

    /**
     * Creates a pretraining domain descriptor
     * @return Pretraining domain descriptor from createpretrainingDomainDescriptor
     */
    @Override
    protected DomainDescriptor<RecordType> domainDescriptor() {
        Properties pretrainingProperties = pretrainingProperties();
        return new PretrainingDomainDescriptor<>(delegate.domainDescriptor(),
                createArguments(), pretrainingFeatureMapperClassName(), pretrainingLabelMapperClassName(),
                pretrainingProperties);
    }

    @Override
    public Properties getReaderProperties(String trainingSet) throws IOException {
        Properties delegateProperties = delegate.getReaderProperties(trainingSet);
        delegateProperties.setProperty("delegate.trainModel", delegate.getClass().getCanonicalName());
        delegateProperties.setProperty("delegate.domainDescriptor",
                delegate.domainDescriptor().getClass().getCanonicalName());
        return delegate.getReaderProperties(trainingSet);
    }

    @Override
    public TrainingArguments createArguments() {
        return delegate.createArguments();
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
     * Properties object with any relevant properties needed to create pretraining feature and label mappers
     * @return properties object
     */
    public abstract Properties pretrainingProperties();
}
