package org.campagnelab.dl.genotype.helpers;

import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.goby.algorithmic.indels.EquivalentIndelRegion;
import org.campagnelab.goby.alignments.processors.ObservedIndel;
import org.campagnelab.goby.algorithmic.algorithm.EquivalentIndelRegionCalculator;


import java.util.Iterator;
import java.util.Set;

/**
 * Created by fac2003 on 12/25/16.
 */
public class GenotypeHelper {


    private EquivalentIndelRegionCalculator equivalentIndelRegionCalculator;



    public static boolean isVariant(String trueGenotype, String referenceBase) {
        return isVariant(true, trueGenotype, referenceBase);
    }


    public static boolean isVariant(Set<String> genotype, String referenceBase) {
        return isVariant(true, genotype, referenceBase);
    }

    public static boolean isNoCall(String genotype) {
        return ("N".equals(genotype) || "N|N".equals(genotype) || "N/N".equals(genotype));
    }



    public static boolean isVariant(boolean considerIndels, String trueGenotype, String reference){
        return isVariant(considerIndels,getAlleles(trueGenotype),reference);
    }

    public static boolean isVariant(boolean considerIndels, Set<String> genotypeSet, String reference) {



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

    public static String fromAlleles(Set<String> alleles){
        StringBuffer sb = new StringBuffer();
        if (alleles.size() > 1){
            for (String allele : alleles){
                sb.append(allele + "|");
            }
            return sb.substring(0,sb.length()-1);
        } else {
            for (String allele : alleles){
                return (allele + "|" + allele).toUpperCase();
            }
            return ".|.";
        }

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

    public static boolean genotypeHasAlleleOrIndel(String trueGenotype, String toSequence, String trueFrom, String fromSequence) {
        boolean hasTo = false;
            Set<String> alleles = getAlleles(trueGenotype);

            for (String allele : getAlleles(trueGenotype)){
                //handle tail bug and append to true genotype (sometimes the goby realignment contains some extra characters at the end
//                if (allele.length() > 1 && toSequence.length() > allele.length() && toSequence.charAt(toSequence.length()-1)!='-') {
//                    //but don't append these extra characters if they are part of an insertion
//                    if (fromSequence.length() >= toSequence.length() && (!fromSequence.substring(allele.length(),fromSequence.length()).contains("-"))){
//                        String altAllele1 = allele + toSequence.substring(allele.length(), toSequence.length()); // append to true genotype
//                        alleles.add(altAllele1);
//                    }
//                }
////                handle additional trailing dash bug where true gentotype realign lacks one trailing dash
//                if (toSequence.charAt(toSequence.length()-1)=='-' && allele.charAt(allele.length()-1)=='-' && toSequence.length() == allele.length()+1){
//                    alleles.add(allele + "-");
//                }


//                if (toSequence.length() > 1 && allele.length() > toSequence.length() && fromSequence.equals(trueFrom.substring(0,fromSequence.length())) && allele.charAt(toSequence.length())!='-'){
//                    String altAllele2 = allele.substring(0,toSequence.length()); //clip 1 to handle flank blug
//                    alleles.add(altAllele2);
//                }
        }
        Iterator<String> iterator = alleles.iterator();
        while (iterator.hasNext()) {
            String oneTrueAllele = iterator.next();
            hasTo |= oneTrueAllele.equals(toSequence);
        }
        return hasTo;
    }


    class callTruePair {
        String trueTo; //true allele
        String ref; //vcf's "from" string
        String from; //goby's from
        String to; //goby's to
        int posRef;
        int referenceIndex;



        /*
        *extend true allele to match length of its ref. eg: ATC A -> ATC A--
        * This is a required step before using observeIndels, which requires placeholder dashes
         */
        void padTrueAndRef(){
            int maxLen = Math.max(trueTo.length(), ref.length());
            trueTo = pad(trueTo,maxLen);
            ref = pad(ref,maxLen);
        }

        void trueToToEquivalent(){

            if (trueTo.length() < 2 || ref.length() < 2){
                //snp encountered, in an indel case one will be longer and the other should have been padded.
                return;
            }
            String trueToAffix = trueTo.substring(1);
            String refAffix = ref.substring(1);


            ObservedIndel indel = new ObservedIndel(posRef, posRef, refAffix, "T");
            EquivalentIndelRegion result = equivalentIndelRegionCalculator.determine(3, indel);



//            assertEquals(3, result.referenceIndex);
//            assertEquals(4, result.startPosition);
//            assertEquals(7, result.endPosition);
//            assertEquals("-TT", result.from);
//            assertEquals("TTT", result.to);
//            assertEquals("AAAC", result.flankLeft);
//            assertEquals("GGGG", result.flankRight);
//            assertEquals("AAAC-TTGGGG", result.fromInContext());
//            assertEquals("AAACTTTGGGG", result.toInContext());
        }








    }

}
