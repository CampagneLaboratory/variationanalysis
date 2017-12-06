package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.performance.SegmentGenotypePredictionTest;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Represents predicted genotypes for a segment.
 */
public class SegmentPrediction extends Prediction {


    SegmentInformationRecords.ReferencePosition startPosition;
    SegmentInformationRecords.ReferencePosition endPosition;
    SegmentGenotypePrediction genotypes;
    IntSet indelPositions = new IntArraySet();

    // predicted colors, one per allele ( up to ploidy), predictedColors[alleleIndex][baseIndex].
    // note that alleleIndex identifies the allele in the predicted genotype. 0 is first allele in 0/1/2,
    // 1 second, and so on.
    //byte[][] predictedColors;

    /**
     * Meta data, either populated from the segment record, or from meta-data in the cache.
     */
    // isVariant[baseIndex] is true when the base does not match the reference.
    boolean[] isVariant;
    // isIndelPosition[baseIndex] is true when the base participates to an indel.
    boolean[] isIndel;

    public SegmentPrediction(SegmentInformationRecords.ReferencePosition startPosition,
                             SegmentInformationRecords.ReferencePosition endPosition,
                             SegmentGenotypePrediction segmentGenotypePrediction) {
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.genotypes = segmentGenotypePrediction;
    }

    public int numPredictionsWhere(SegmentGenotypePredictionTest isTrueForBaseIndex) {
        int count = 0;
        final int numBases = numBases();
        for (int baseIndex = 0; baseIndex < numBases; baseIndex++) {
            if (isTrueForBaseIndex.test(this.genotypes,baseIndex)) {
                count+=1;
            }
        }

        return count;
    }

    public int numBases() {
        return genotypes.numBases();
    }

    public int getStartPosition() {
        return startPosition.getLocation();
    }

    public int getEndPosition() {
        return endPosition.getLocation();
    }
    
    public SegmentGenotypePrediction getGenotypes() {
        return genotypes;
    }

    public void inspectRecord(SegmentInformationRecords.SegmentInformation record) {
        int baseIndex = 0;
        for (SegmentInformationRecords.Sample sample : record.getSampleList()) {
            int previousLocation = 0;
            for (SegmentInformationRecords.Base base :sample.getBaseList()) {
                if (base.getLocation() == previousLocation || genotypes.predictedGenotypes[baseIndex].contains("-") ) {
                    //this is an indel (at least, two bases at the same location with a predicted gap)
                     this.indelPositions.add(base.getLocation());
                }
                previousLocation = base.getLocation();
                baseIndex++;
            }
        }


    }

    public Set<String> predictedAllelesForIndels(int position, int howMany) {
        Set<String> alleles = new HashSet<>();
        String firstAllele = "";
        String secondAllele = "";
        int i = 0;
        while (i < howMany && i< genotypes.predictedGenotypes.length && genotypes.predictedGenotypes[i] != null) {
            String[] parts = genotypes.predictedGenotypes[position + i++].split("");
            firstAllele += parts[0];
            secondAllele += parts[1];
        }
        alleles.add(firstAllele);
        alleles.add(secondAllele);
        return alleles;
    }

    public Set<String> predictedAllelesForSnps(int position) {
        Set<String> alleles = new HashSet<>();
        String[] genotypesAtBase = genotypes.predictedGenotypes[position].split("");
        alleles.addAll(Arrays.asList(genotypesAtBase));
        return alleles;
    }


    public Set<String> trueAlleles(int position) {
        Set<String> alleles = new HashSet<>();
        alleles.addAll(Arrays.asList(genotypes.trueGenotypes[position].split("")));
        return alleles;
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

    /**
     * Checks whether the base is at a predicted indel positon.
     * @param base
     * @return
     */
    public boolean isIndelPosition(SegmentInformationRecords.Base base) {
        return this.indelPositions.contains(base.getLocation());
    }

    public String referenceAllelesForIndels(ObjectList<SegmentInformationRecords.Base> basesAt) {
        final String[] reference = {""};
        basesAt.forEach(base -> reference[0] += base.getReferenceAllele());
        return reference[0];
    }
}
