package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictWithVecFile;
import org.campagnelab.dl.genotype.mappers.MetaDataLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.IOException;
import java.util.List;

public class PredictGWithVecFile extends PredictWithVecFile<BaseInformationRecords.BaseInformation> {
    public PredictGWithVecFile(DomainDescriptor domainDescriptor, String vecPath) throws IOException {
        super(domainDescriptor, vecPath);
    }



}
