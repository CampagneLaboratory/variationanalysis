package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.ExportTensors;
import org.campagnelab.dl.somatic.learning.SomaticTrainingArguments;
import org.campagnelab.dl.somatic.learning.domains.SomaticMutationDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

public class ExportTensorsS extends ExportTensors<BaseInformationRecords.BaseInformation> {
    private SomaticMutationDomainDescriptor domainDescriptor;

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor(String featureMapperClassName,
                                                                                        List<String> trainingSets,
                                                                                        int genomicContextLength,
                                                                                        float labelSmoothingEpsilon,
                                                                                        int ploidy,
                                                                                        int extraGenotypes) {
        SomaticTrainingArguments trainingArguments=new SomaticTrainingArguments();
        trainingArguments.featureMapperClassname=featureMapperClassName;

        trainingArguments.genomicContextLength=genomicContextLength;
        trainingArguments.labelSmoothingEpsilon=labelSmoothingEpsilon;
        trainingArguments.ploidy=ploidy;
        trainingArguments.extraGenotypes=extraGenotypes;
        trainingArguments.trainingSets=new ArrayList<>();
        trainingArguments.trainingSets.addAll(trainingSets);

        SomaticMutationDomainDescriptor domainDescriptor=new SomaticMutationDomainDescriptor(trainingArguments);
        this.domainDescriptor = domainDescriptor;
        return domainDescriptor;
    }



    /**
     * Record arguments to the properties, that need to be provided to feature/label mappers.
     *
     * @param properties
     */
    @Override
    protected void decorateProperties(Properties properties) {
        if (domainDescriptor != null) {
            domainDescriptor.decorateProperties(properties);
        }
    }

    @Override
    public Properties getReaderProperties(String trainingSet) throws IOException {
        try (SequenceBaseInformationReader reader = new SequenceBaseInformationReader(trainingSet)) {
            final Properties properties = reader.getProperties();
            reader.close();
            return properties;
        }
    }

    public static void main(String[] args) {

        ExportTensorsS exportTensorsS = new ExportTensorsS();
        exportTensorsS.parseArguments(args, "ExportTensorS", exportTensorsS.createArguments());
        exportTensorsS.execute();
    }
}
