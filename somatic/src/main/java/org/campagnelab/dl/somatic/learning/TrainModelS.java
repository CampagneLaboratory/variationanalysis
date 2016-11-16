package org.campagnelab.dl.somatic.learning;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.TrainModel;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.somatic.learning.domains.SomaticMutationDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Properties;

/**
 * Train Somatic implemented with the Generic TrainModel
 */
public class TrainModelS extends TrainModel<BaseInformationRecords.BaseInformation> {
    static private Logger LOG = LoggerFactory.getLogger(TrainModelS.class);

    public static void main(String[] args) {

        TrainModelS tool = new TrainModelS();
        tool.parseArguments(args, "TrainModelS", tool.createArguments());
        if (tool.args().trainingSets.size() == 0) {
            System.out.println("Please add exactly one training set to the args().");
            return;
        }
        assert !tool.args().errorEnrichment : "This tool does not support error enrichment";
        tool.execute();
        tool.writeModelingConditions(tool.getRecordingArguments());
    }

    @Override
    public TrainingArguments createArguments() {
        return new SomaticTrainingArguments();
    }

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor() {
        return new SomaticMutationDomainDescriptor((SomaticTrainingArguments) args());
    }

    @Override
    public Properties getReaderProperties(String trainingSet) throws IOException {

        SequenceBaseInformationReader reader = new SequenceBaseInformationReader(trainingSet);
        final Properties properties = reader.getProperties();
        reader.close();
        return properties;

    }
}
