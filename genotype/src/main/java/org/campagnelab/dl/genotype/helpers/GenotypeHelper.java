package org.campagnelab.dl.genotype.helpers;

import org.campagnelab.dl.genotype.predictions.GenotypePrediction;

import java.util.Set;

/**
 * Created by fac2003 on 12/25/16.
 */
public class GenotypeHelper {
    public static boolean isVariant(String genotype, String referenceBase) {
        return isVariant(true, genotype, referenceBase);
    }

    public static boolean isVariant(boolean considerIndels, String genotype, String reference) {

        Set<String> genotypeSet = getAlleles(genotype);
        //handle special case where a homozygous deletion like AG->A is not caught and matches ref
        //such cases chould not be counted as variants.
        boolean matchesRef = false;

        if (genotypeSet.size() == 1 && genotypeSet.iterator().next().equals(reference)) {
            matchesRef = true;
        }
        if (genotypeSet.size() > 1) {
            if (genotypeSet.stream().map(allele -> allele.length()).distinct().count() == 1) {
                // alleles have same length
                if (considerIndels) {
                    return true;
                } else {
                    String referenceBase=reference.substring(0,1);
                 return!   referenceBase.equals(genotypeSet.stream().map(allele -> Character.toString(allele.charAt(0))).distinct().findFirst().get());
                }
            }
            matchesRef = false;
        }
        return !matchesRef;
    }

    public static Set<String> getAlleles(String genotype) {
        String trueGenotype = genotype.toUpperCase();
    return    GenotypePrediction.alleles(trueGenotype);
    }

    public static boolean isIndel(String genotype) {
        if (genotype != null) {
            return (genotype.contains("-"));
        }
        return false;
    }
}