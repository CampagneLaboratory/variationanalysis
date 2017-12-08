package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.performance.SegmentGenotypePredictionTest;
import org.campagnelab.dl.genotype.tools.VCFLine;
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
    IntSet hasGaps = new IntArraySet();

    // predicted colors, one per allele ( up to ploidy), predictedColors[alleleIndex][baseIndex].
    // note that alleleIndex identifies the allele in the predicted genotype. 0 is first allele in 0/1/2,
    // 1 second, and so on.
    //byte[][] predictedColors;

    /**
     * Meta data, either populated from the segment record, or from meta-data in the cache.
     */
    // isVariant[baseIndex] is true when the base does not match the reference.
    boolean[] isVariant;

    private String referenceId;

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
        referenceId = record.getStartPosition().getReferenceId();
        int baseIndex = 0;
        for (SegmentInformationRecords.Sample sample : record.getSampleList()) {
            int previousLocation = 0;
            for (SegmentInformationRecords.Base base :sample.getBaseList()) {
                if (base.getLocation() == previousLocation || genotypes.predictedGenotypes[baseIndex].contains("-") ) {
                    //this is an indel (at least, two bases at the same location with a predicted gap)
                     this.hasGaps.add(base.getLocation());
                }
                previousLocation = base.getLocation();
                baseIndex++;
            }
        }


    }

    public boolean hasPredictedGap(int segmentIndex) {
      return genotypes.predictedGenotypes[segmentIndex].contains("-");
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


    public String referenceAlleles(VCFLine line) {
        final String[] reference = {""};
        line.forEach(base -> reference[0] += base.getKey().getReferenceAllele());
        return reference[0];
    }

    public Set<String> predictedAlleles(VCFLine line) {
        Set<String> alleles = new HashSet<>();
        String firstAllele = "";
        String secondAllele = "";
        for (VCFLine.IndexedBase base : line) {
            String[] parts = genotypes.predictedGenotypes[base.getValue()].split("");
            firstAllele += parts[0];
            secondAllele += parts[1];
        }
        alleles.add(firstAllele);
        alleles.add(secondAllele);
        return alleles;
    }


    public String getReferenceId() {
        return referenceId;
    }
}
