package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.ExportTensors;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.ArrayList;
import java.util.List;

public class ExportTensorsG extends ExportTensors<BaseInformationRecords.BaseInformation> {
    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor(String featureMapperClassName,
                                                                                        List<String> trainingSets) {
        GenotypeTrainingArguments trainingArguments=new GenotypeTrainingArguments();
        trainingArguments.featureMapperClassname=featureMapperClassName;
        trainingArguments.trainingSets=new ArrayList<>();

        trainingArguments.trainingSets.addAll(trainingSets);

        GenotypeDomainDescriptor domainDescriptor=new GenotypeDomainDescriptor(trainingArguments);
        return domainDescriptor;
    }

    public static void main(String[] args) {

        ExportTensorsG exportG = new ExportTensorsG();
        exportG.parseArguments(args, "ExportTensorG", exportG.createArguments());
        exportG.execute();
    }
}
