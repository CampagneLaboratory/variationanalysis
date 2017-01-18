package org.campagnelab.dl.genotype.helpers;

import org.campagnelab.dl.genotype.predictions.GenotypePrediction;

import java.util.Iterator;
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
            if (considerIndels){
                matchesRef = (reference.equals(allele));
            } else {
                //just check first base in this case.
                if (reference.charAt(0)==(allele.charAt(0))) {
                    matchesRef = true;
                }
            }

        }
        if (genotypeSet.size() > 1) {
            if (considerIndels) {
                return true;
            } else {
                String referenceBase = reference.substring(0, 1);
                return !referenceBase.equals(genotypeSet.stream().map(allele -> Character.toString(allele.charAt(0))).distinct().findFirst().get());
            }
        }
        return !matchesRef;
    }

    public static Set<String> getAlleles(String genotype) {
        String trueGenotype = genotype.toUpperCase();
        return GenotypePrediction.alleles(trueGenotype);
    }


    /**
     * genotype must have two alleles
     * @param reference
     * @param genotype
     * @return
     */
    public static boolean isIndel(String reference, String genotype) {
        if (genotype != null) {
            return (genotype.length()>3 || reference.length()>1);
        }
        return false;
    }

    public static boolean matchingGenotypesWithN(String a, String b) {
        if (isNoCall(a) || isNoCall(b)) {
            return true;
        }
      return  matchingGenotypes(a,b);
    }

    public static boolean matchingGenotypes(String a, String b) {

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

    /**
     * Return true iff the true genotype has an allele matching toSequence.
     * @param trueGenotype a true genotype, e.g., A/T
     * @param testAllele the sequence of an allele, e.g., A
     * @return True if and only if one of the true genotype alleles is matching testAllele.
     */
    public static boolean genotypeHasAllele(String trueGenotype, String testAllele) {
        Set<String> alleles = getAlleles(trueGenotype);
        Iterator<String> iterator = alleles.iterator();
        while (iterator.hasNext()) {
            String oneTrueAllele = iterator.next();
            if (matchingGenotypes(oneTrueAllele, testAllele)) {
                return true;
            }
        }
        return false;
    }



    public static String pad(String s, int maxLength){
        StringBuffer toPad = new StringBuffer(s);
        for (int i = 1; i <= maxLength; i++) {
            if (toPad.length() < maxLength) {
                toPad.append("-");
            }
        }
        return toPad.toString();
    }

    public static int maxLength (Set<String> strings){
        int max = 0;
        for (String s : strings){
            if (s.length() > max) {
                max = s.length();
            }
        }
        return max;
    }


    public static String padMulti(String genotype, int maxLength){
        StringBuffer toPad = new StringBuffer();
        Set<String> alleles = getAlleles(genotype);
        for (String allele : alleles) {
            toPad.append(pad(allele, maxLength) + "|");
        }
        //remove final |
        toPad.deleteCharAt(toPad.length()-1);
        return toPad.toString();
    }

    public static boolean genotypeHasIndel(String trueGenotype, String toSequence, String trueFrom, String fromSequence) {
        boolean hasTo = genotypeHasAllele(trueGenotype,toSequence);
        boolean fromMatches = matchingGenotypes(trueFrom,fromSequence);
        return hasTo && fromMatches;
    }
}
