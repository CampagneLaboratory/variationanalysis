package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.ExportTensors;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

public class ExportTensorsG extends ExportTensors<BaseInformationRecords.BaseInformation> {
    private GenotypeDomainDescriptor genotpeDomainDescriptor;

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor(String featureMapperClassName,
                                                                                        List<String> trainingSets) {
        GenotypeTrainingArguments trainingArguments=new GenotypeTrainingArguments();
        trainingArguments.featureMapperClassname=featureMapperClassName;
        trainingArguments.trainingSets=new ArrayList<>();

        trainingArguments.trainingSets.addAll(trainingSets);

        GenotypeDomainDescriptor domainDescriptor=new GenotypeDomainDescriptor(trainingArguments);
        this.genotpeDomainDescriptor = domainDescriptor;
        return domainDescriptor;
    }

    /**
     * Record arguments to the properties, that need to be provided to feature/label mappers.
     *
     * @param properties
     */
    @Override
    protected void decorateProperties(Properties properties) {
        if (genotpeDomainDescriptor != null) {
            genotpeDomainDescriptor.decorateProperties(properties);
        }
    }

    public static void main(String[] args) {

        ExportTensorsG exportG = new ExportTensorsG();
        exportG.parseArguments(args, "ExportTensorG", exportG.createArguments());
        exportG.execute();
    }
}
