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
        return new PretrainingDomainDescriptor<>(delegate.domainDescriptor(), recordToSequenceLength(),
                createArguments());
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
     * Create a function that takes in a sequence record and returns its length
     * Should preferably be an O(1) function
     * @return Function to convert sequence record to length
     */
    public abstract Function<RecordType, Integer> recordToSequenceLength();
}
