package org.campagnelab.dl.genotype.segments;


import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.genotype.mappers.SingleBaseLabelMapperV1;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;


/**
 * Created by mas2182 on 10/17/17.
 */
public class SegmentLabelMapper {

    private final int ploidy;
    private final int numberOfLabelsPerBase; // Combinations of  A,C,T,G,-
    private final Set<String> labels;
    private final char[] alleles = new char[]{'A','C','T','G','-'};

    static private Logger LOG = LoggerFactory.getLogger(SegmentLabelMapper.class);

    public static void main(String[] args) {
        new SegmentLabelMapper(3);
    }

    public SegmentLabelMapper(int ploidy) {
        this.ploidy = ploidy;
        numberOfLabelsPerBase = SingleBaseLabelMapperV1.getNumValues(ploidy, alleles.length);
        labels = new ObjectArraySet<>(numberOfLabelsPerBase);
        this.labels.addAll(this.buildLabels(ploidy));
        this.buildLabelMap();
        //IndexedIdentifier elementLabels = new IndexedIdentifier(numberOfLabelsPerBase);
       // final MutableString elementLabel = new MutableString(label).compact();
       // final int elementIndex = elementLabels.registerIdentifier(elementLabel);

    }

    private void buildLabelMap() {

    }

    /**
     * Creates all the possible combinations with the alleles.
     */
    private List<String> buildLabels(int deep) {
        List<String> combinedLabels = new ObjectArrayList<>();
            for (int from = 0; from < alleles.length; from++) {
                if (deep == 2) {
                    for (int index = from; index < alleles.length ; index++) {
                        combinedLabels.add(sortAlphabetically(String.format("%c%c", alleles[from], alleles[index])));
                    }
                }  else {
                    int finalFrom = from;
                    buildLabels(deep-1).forEach(label -> {
                        combinedLabels.add(sortAlphabetically(String.format("%c%s",alleles[finalFrom],label)));
                    });

                }
            }

        return combinedLabels;
    }

    private String sortAlphabetically(String toSort) {
        char[] elements = toSort.toCharArray();
        Arrays.sort(elements);
        return new String(elements);
    }
    
    public int numberOfLabels() {
        return numberOfLabelsPerBase;
    }

}
