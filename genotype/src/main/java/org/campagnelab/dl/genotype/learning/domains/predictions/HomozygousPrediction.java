package org.campagnelab.dl.genotype.learning.domains.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Set;

/**
 * Created by rct66 on 11/12/16.
 */
public class HomozygousPrediction extends Prediction {

    public Set<String> trueGenotype;
    public String trueGenotypeFormat;
    public String predictedHomozygousGenotype;
    public double probability;

    public <BaseInformation> void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {
        trueGenotype = getGenotype(currentRecord);
        trueGenotypeFormat = (trueGenotype.size()==0)?"./.":currentRecord.getTrueGenotype();
        return;
    }

    public Set<String> getGenotype(BaseInformationRecords.BaseInformation currentRecord){
        Set<String> alleles = new ObjectArraySet<>();
        for (String allele : currentRecord.getTrueGenotype().split("/")){
            alleles.add(allele);
        }
        alleles.remove("");
        alleles.remove("?");
        alleles.remove(".");
        return alleles;
    }
}