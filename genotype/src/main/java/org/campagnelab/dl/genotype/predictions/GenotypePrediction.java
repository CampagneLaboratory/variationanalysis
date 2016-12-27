package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.Set;

/**
 * Describes a genotype prediction. Helper method set aggregates individual model predictions.
 * Created by fac2003 on 12/18/16.
 */
public class GenotypePrediction extends Prediction {

    /**
     * Genotype called by the model.
     */
    public String predictedGenotype;
    /**
     * Genotype called by the model.
     */
    public String trueGenotype;
    /**
     * The probability of the called genotype according to the model. Forecast probability.
     */
    public double overallProbability;

    /**
     * Whether the true genotypes contains at least one non-reference allele.
     */
    public boolean isVariant;
    /**
     * Whether the true genotypes contains at least one indel.
     */
    public boolean isIndel;

    /**
     * Indicates the confidence that the genotype is a variant.
     */
    public double isVariantProbability;

    public GenotypePrediction(String predictedGenotype, String trueGenotype) {
        this.predictedGenotype = predictedGenotype;
        this.trueGenotype = trueGenotype;
    }

    public GenotypePrediction() {
    }

    /**
     * Inspect a record to obtain the true genotype and/or additional meta-data about the record.
     */
    public void inspectRecord(BaseInformationRecords.BaseInformation currentRecord) {
        trueGenotype = currentRecord.getTrueGenotype();
        isVariant = currentRecord.getSamples(0).getIsVariant();
        isIndel = trueGenotype.contains("-");
    }


    /**
     * Returns whether or not the prediction is correct, will be used for statistics.
     */
    public boolean isCorrect() {
        return GenotypeHelper.matchingGenotypes(trueGenotype, predictedGenotype);
    }

    public boolean isVariant() {
        return isVariant;
    }

    public boolean isIndel() {
        return isIndel;
    }

    public Set<String> predictedAlleles() {
        return alleles(predictedGenotype);
    }

    public Set<String> trueAlleles() {
        return alleles(trueGenotype);
    }

    public static Set<String> alleles(String genotype) {
        genotype = genotype.toUpperCase();
        ObjectSet<String> result = new ObjectArraySet<>();
        Collections.addAll(result, genotype.split("[|/]"));
        result.remove("|");
        result.remove("/");
        result.remove("?");
        result.remove(".");
        result.remove("");
        return result;
    }

}
