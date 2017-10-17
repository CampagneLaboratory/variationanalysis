package org.campagnelab.dl.genotype.segments;


import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.lang.MutableString;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;


/**
 * Created by mas2182 on 10/17/17.
 */
public class SegmentLabelMapper {

    private final int ploidy;
    private final int numberOfLabelsPerBase; // Combinations of  A,C,T,G,-
    private final Set<String> labels;
    private final char[] alleles = new char[]{'A','C','T','G','-'};
    private final IndexedIdentifier indexedLabels;
    static private Logger LOG = LoggerFactory.getLogger(SegmentLabelMapper.class);

    public static void main(String[] args) {
        SegmentLabelMapper mapper = new SegmentLabelMapper(3);
        Properties props = new Properties();
        mapper.writeMap(props);
        try(OutputStream output = new FileOutputStream("indexedLabels.properties")) {
            props.store(output, null);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public SegmentLabelMapper(int ploidy) {
        this.ploidy = ploidy;
        this.labels = new ObjectArraySet<>();
        this.labels.addAll(this.buildLabels(ploidy));
        this.numberOfLabelsPerBase = this.labels.size();
        this.indexedLabels = this.buildLabelMap();
    }

    /**
     * Creates a map <index, label>
     * @return
     */
    private IndexedIdentifier buildLabelMap() {
        IndexedIdentifier indexedLabels = new IndexedIdentifier(numberOfLabelsPerBase);
        for (String label : this.labels) {
            final MutableString elementLabel = new MutableString(label).compact();
            final int elementIndex = indexedLabels.registerIdentifier(elementLabel);
            System.out.println(String.format("%d -> %s", elementIndex, label));
        }
        return indexedLabels;
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

    public void writeMap(final Properties properties) {
        for (Map.Entry<MutableString, Integer> entry : this.indexedLabels.entrySet()) {
            properties.put(Integer.toString(entry.getValue()), entry.getKey().toString());
        }
    }

    /**
     * Sorts the chars in the list in alphabetic order.
     * @param toSort
     * @return
     */
    private String sortAlphabetically(String toSort) {
        char[] elements = toSort.toCharArray();
        Arrays.sort(elements);
        return new String(elements);
    }
    
    public int numberOfLabels() {
        return numberOfLabelsPerBase;
    }

}
