package org.campagnelab.dl.genotype.learning.domains.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.predictions.AbstractGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Set;

/**
 * Created by rct66 on 11/12/16.
 * Homozygous refers only prediction, not the true label.
 */
public class CombinedPrediction extends AbstractGenotypePrediction {

    //todo: fill these bools, they will be used for statistics like sensitivity.
    public boolean isVariant;
    public boolean isIndel;
    public Set<String> trueGenotype;
    public String trueGenotypeFormat;
    public double probability;
    /**
     * True when the model predicts the site to be homozygous. False when the model predicts het.
     */
    public void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {
        trueGenotype = getGenotype(currentRecord);
        trueGenotypeFormat = (trueGenotype.size()==0)?"./.":currentRecord.getTrueGenotype();
        assert trueGenotypeFormat!="./.":"record with no true genotype label encountered";
        isVariant = currentRecord.getSamples(0).hasIsVariant() && currentRecord.getSamples(0).getIsVariant();
        isIndel = trueGenotype.contains("-");
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

    public boolean isCorrect(){
        return alleles(trueGenotypeFormat.toUpperCase()).equals(alleles(predictedGenotype.toUpperCase()));
    }
}