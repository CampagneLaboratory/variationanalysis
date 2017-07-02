package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.algorithmic.dsv.SampleCountInfo;
import org.campagnelab.goby.predictions.MergeIndelFrom;
import org.campagnelab.goby.util.Variant;

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
     * The 'from' field corresponding to predicted call
     */
    public String predictedFrom;
    /**
     * The 'from' field corresponding to true call
     */
    public String trueFrom;
    /**
     * The index of the reference allele in the goby counts (sbi), or -1 if no count matched the reference.
     */
    public int referenceGobyIndex;
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
     * Whether the predicted genotypes contain at least one indel.
     */
    public boolean isPredictedIndel;

    /**
     * Indicates the confidence that the genotype is a variant.
     */
    public double isVariantProbability;
    private boolean predictedSNP;

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

        isVariant = currentRecord.getSamples(0).getIsVariant();
        Set<Variant.FromTo> predictedFromTos = new ObjectArraySet<>(2);
        Set<Variant.FromTo> trueFromTos = new ObjectArraySet<>(2);
        for (BaseInformationRecords.CountInfo c : currentRecord.getSamples(0).getCountsList()) {
            if (predictedAlleles().contains(c.getToSequence())) {
                predictedFromTos.add(new Variant.FromTo(c.getFromSequence(), c.getToSequence()));
            }
            if (c.getIsCalled()) {
                trueFromTos.add(new Variant.FromTo(c.getFromSequence(), c.getToSequence()));
            }
        }
        MergeIndelFrom mergedTrue = new MergeIndelFrom(trueFromTos);
        MergeIndelFrom mergedPrediction = new MergeIndelFrom(predictedFromTos);
        trueGenotype = GenotypeHelper.fromAlleles(mergedTrue.getTos());
        predictedGenotype = GenotypeHelper.fromAlleles(mergedPrediction.getTos());
        trueFrom = mergedTrue.getFrom();
        if (predictedGenotype.equals(".|.")) {
            predictedFrom = currentRecord.getReferenceBase();
        } else {
            predictedFrom = mergedPrediction.getFrom();
        }
        isIndel = GenotypeHelper.isIndel(mergedTrue);
        isPredictedIndel = GenotypeHelper.isIndel(mergedPrediction);

        /*
       //we need to check from and to fields for a genotype greater than length 1
        //can't use existence of a dash "-" because some indels don't use them
        isIndel = false;
        for (String trueTo : trueAlleles()) {
            isIndel |= trueTo.length() > 1;
        }
        isIndel |= trueFrom.length() > 1;

        isPredictedIndel = false;
        for (String predTo : predictedAlleles()) {
            isPredictedIndel |= predTo.length() > 1;
            isPredictedIndel |= predTo.contains("-");
        }
        isPredictedIndel |= predictedFrom.length() > 1;
        isPredictedIndel |= predictedFrom.contains("-");
        */
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

    public boolean isPredictedIndel() {
        return isPredictedIndel;
    }

    public boolean isSnp() {
        return !isIndel && isVariant;
    }

    public boolean isPredictedSnp() {
        return !isPredictedIndel;
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



    public void set(MetadataPrediction metaData) {
        this.isIndel = metaData.isIndel;
        this.isVariant = metaData.isVariant;
        this.referenceGobyIndex = metaData.referenceGobyIndex;
    }

    public static String DECISION_THRESHOLD_PROPERTY = "genotypes.decisionThreshold";
}
