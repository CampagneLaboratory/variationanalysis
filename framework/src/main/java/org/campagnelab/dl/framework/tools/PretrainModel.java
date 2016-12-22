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
    private DomainDescriptor<RecordType> delegateDescriptor;

    /**
     * Create a pretraining model by using another TrainModel as a delegate
     * @param delegate delegate TrainModel
     */
    public PretrainModel(TrainModel<RecordType> delegate) {
        this.delegate = delegate;
        delegateDescriptor = delegate.domainDescriptor();
    }


    @Override
    protected DomainDescriptor<RecordType> domainDescriptor() {
        return createPretrainingDomainDescriptor();
    }

    @Override
    public Properties getReaderProperties(String trainingSet) throws IOException {
        Properties delegateProperties = delegate.getReaderProperties(trainingSet);
        delegateProperties.setProperty("delegate.train_Model", delegate.getClass().getCanonicalName());
        delegateProperties.setProperty("delegate.domain_descriptor",
                delegate.domainDescriptor().getClass().getCanonicalName());
        String eosString = args().eosIndex == null ? "null" : Integer.toString(args().eosIndex);
        delegateProperties.setProperty("delegate.eos_index", eosString);
        return delegateProperties;
    }

    @Override
    public TrainingArguments createArguments() {
        return delegate.createArguments();
    }

    /**
     * Get the delegate TrainModel
     * @return delegate TrainModel
     */
    public TrainModel<RecordType> delegateTrainModel() {
        return delegate;
    }

    /**
     * Get the domain of the delegate TrainModel
     * @return domain of delegate TrainModel
     */
    public DomainDescriptor<RecordType> delegateDomainDescriptor() {
        return delegateDescriptor;
    }

    /**
     * Create a pretraining domain descriptor to be used as PretrainModel's domain descriptor
     * @return Pretraining domain descriptor
     */
    public abstract PretrainingDomainDescriptor<RecordType> createPretrainingDomainDescriptor();
}
