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

    public static boolean isNoCall(String genotype) {
        return ("N".equals(genotype) || "N|N".equals(genotype) || "N/N".equals(genotype));
    }

    public static boolean isVariant(boolean considerIndels, String genotype, String reference) {

        Set<String> genotypeSet = getAlleles(genotype);
        //handle special case where a homozygous deletion like AG->A is not caught and matches ref
        //such cases chould not be counted as variants.
        boolean matchesRef = false;


        if (genotypeSet.size() == 1) {
            String allele = genotypeSet.iterator().next();
            // When there is only one allele, only the first bases need to match with the reference,
            // up to the length of the reference:
            // i.e., genotype: TC/TC ref: T, or genotype: TCA/TCA ref: TC
            int refLength = reference.length();
            if (reference.substring(0, refLength).equals(allele.substring(0, refLength))) {
                matchesRef = true;
            }
        }
        if (genotypeSet.size() > 1) {
            if (genotypeSet.stream().map(allele -> allele.length()).distinct().count() == 1) {
                // alleles have same length
                if (considerIndels) {
                    return true;
                } else {
                    String referenceBase = reference.substring(0, 1);
                    return !referenceBase.equals(genotypeSet.stream().map(allele -> Character.toString(allele.charAt(0))).distinct().findFirst().get());
                }
            }
            matchesRef = false;
        }
        return !matchesRef;
    }

    public static Set<String> getAlleles(String genotype) {
        String trueGenotype = genotype.toUpperCase();
        return GenotypePrediction.alleles(trueGenotype);
    }

    public static boolean isIndel(String genotype) {
        if (genotype != null) {
            return (genotype.contains("-"));
        }
        return false;
    }

    public static boolean matchingGenotypes(String a, String b) {
        if (isNoCall(a) || isNoCall(b)) {
            return true;
        }
        Set<String> allelesA = getAlleles(a);
        Set<String> allelesB = getAlleles(b);
        boolean allelesMatch = allelesA.equals(allelesB);
        if (allelesMatch) {
            return true;
        }
        if (allelesA.size() == 1 && allelesB.size() == 1) {

            String oneA = allelesA.iterator().next();
            String oneB = allelesB.iterator().next();
            if (a.contains("-") || b.contains("-")) {
             // do not allow prefix match for indels.
                return false;
            }
            if (oneA.length() > oneB.length()) {
                return oneA.startsWith(oneB);
            }
            if (oneB.length() > oneA.length()) {
                return oneB.startsWith(oneA);
            }
        }
        return false;
    }
}
