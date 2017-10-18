package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.TrainModel;
import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.genotype.learning.domains.GenotypeSegmentDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Properties;

/**
 * Train a Genotype model, from .ssi files. Implemented by specializing the framework TrainModel class.
 */
public class TrainModelGS extends TrainModel<SegmentInformationRecords.SegmentInformation> {
    static private Logger LOG = LoggerFactory.getLogger(TrainModelGS.class);

    public static void main(String[] args) {

        TrainModelGS tool = new TrainModelGS();
        tool.parseArguments(args, "TrainModelGS", tool.createArguments());
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
        return new SegmentTrainingArguments();
    }

    @Override
    protected DomainDescriptor<SegmentInformationRecords.SegmentInformation> domainDescriptor() {
        return new GenotypeSegmentDomainDescriptor((SegmentTrainingArguments) args());
    }

    @Override
    public Properties getReaderProperties(String trainingSet) throws IOException {

        return GenotypeSegmentDomainDescriptor.getReaderProperties(trainingSet);
    }


}
