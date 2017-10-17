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
    private final List<String> labels;
    private final char[] alleles = new char[]{'A','C','T','G','-'};

    static private Logger LOG = LoggerFactory.getLogger(SegmentLabelMapper.class);

    public static void main(String[] args) {
        new SegmentLabelMapper(2);
    }

    public SegmentLabelMapper(int ploidy) {
        this.ploidy = ploidy;
        numberOfLabelsPerBase = SingleBaseLabelMapperV1.getNumValues(ploidy, 5);
        labels = new ObjectArrayList<>(numberOfLabelsPerBase);
        this.labels.addAll(this.buildLabels(ploidy));
        Collections.sort(this.labels);
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
        List<String> partialLabels = new ObjectArrayList<>();
            for (int from = 0; from < alleles.length; from++) {
                if (deep == 2) {
                    for (int index = from; index < alleles.length; index++) {
                        partialLabels.add(String.format("%c%c", alleles[from], alleles[index]));
                    }
                }  else {
                    List<String> deepestLabels = buildLabels(deep-1);
                    int finalFrom = from;
                    buildLabels(deep-1).forEach(label -> {
                        partialLabels.add(String.format("%c%s",alleles[finalFrom],label));
                    });

                }
            }

        return partialLabels;
    }

    
    public int numberOfLabels() {
        return numberOfLabelsPerBase;
    }

}
