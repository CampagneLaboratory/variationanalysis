package org.campagnelab.dl.genotype.segments;


import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.lang.MutableString;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * Created by mas2182 on 10/17/17.
 */
public class SegmentLabelMapper {

    private final int ploidy;
    private final int numberOfLabelsPerBase; // Combinations of  A,C,T,G,-
    private final Set<String> labels;
    private final char[] alleles = new char[]{'A', 'C', 'T', 'G', '-', 'N'};
    private final IndexedIdentifier indexedLabels;
    static private Logger LOG = LoggerFactory.getLogger(SegmentLabelMapper.class);

    public static void main(String[] args) {
        SegmentLabelMapper mapper = new SegmentLabelMapper(Integer.valueOf(args[0]));
        float[] position = mapper.map(args[1], Arrays.stream(args[2].split(","))
                .mapToInt(Integer::valueOf).boxed().collect(Collectors.toList()));
        System.out.println(Arrays.toString(position));
        Properties props = new Properties();
        mapper.writeMap(props);
        try (OutputStream output = new FileOutputStream("indexedLabels.properties")) {
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
     *
     * @return
     */
    private IndexedIdentifier buildLabelMap() {
        IndexedIdentifier indexedLabels = new IndexedIdentifier(numberOfLabelsPerBase);
        for (String label : this.labels) {
            final MutableString elementLabel = new MutableString(label).compact();
            final int elementIndex = indexedLabels.registerIdentifier(elementLabel);
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
                for (char to : alleles) {
                    combinedLabels.add(String.format("%c%c", alleles[from], to));
                }
            } else {
                int finalFrom = from;
                buildLabels(deep - 1).forEach(label -> {
                    combinedLabels.add(String.format("%c%s", alleles[finalFrom], label));
                });

            }
        }
        return combinedLabels;
    }

    /**
     * Maps the given chromosome sequence to one hot encoding.
     *
     * @param alleles set of alleles in the format A/B/C where A,B or C are one of A,C,T,G,-
     * @param indices list of gobyTrueGenotypeIndices corresponding to each allele
     * @return a float array of zeros, except for the position of the sequence (which is equals to 1).
     */
    public float[] map(String alleles, List<Integer> indices) {
        float[] position = new float[this.numberOfLabelsPerBase];
        String clean = alleles.replace("/", "");
        if (clean.length() != ploidy)
            throw new IllegalArgumentException(alleles + " is not of the expected length (" + this.ploidy + ")");
        final MutableString toFind;
        // If sizes unequal, generate indices as range 0, 1,...,n+1, n, of same length as alleles
        if (indices.size() != clean.length()) {
            indices = IntStream.range(0, clean.length()).boxed().collect(Collectors.toList());
        }
        toFind = new MutableString(sortByIndices(clean, indices)).compact();
        int foundAt = this.indexedLabels.getInt(toFind);
        assert foundAt >= 0 : String.format("genotype %s not found in map", toFind);
        if (foundAt < position.length - 1 && foundAt >= 0) {
            position[foundAt] = 1;
        } else {
            position[foundAt] = -1;
        }
        return position;
    }

    /**
     * Add the indexes and labels to the properties.
     *
     * @param properties
     */
    public void writeMap(final Properties properties) {
        for (Map.Entry<MutableString, Integer> entry : this.indexedLabels.entrySet()) {
            properties.put("genotype.segment.label." + Integer.toString(entry.getValue()), entry.getKey().toString());
        }
        properties.put("genotype.segment.label.numOfEntries", Integer.toString(this.indexedLabels.entrySet().size()));
    }

    /**
     * Sorts the chars in the list in alphabetic order.
     *
     * @param toSort
     * @return
     */
    private String sortAlphabetically(final String toSort) {
        char[] elements = toSort.toCharArray();
        Arrays.sort(elements);
        return new String(elements);
    }

    static class CharWithIndex {
        private final String allele;
        private final int index;

        CharWithIndex(char allele, int index) {
            this.allele = "" + allele;
            this.index = index;
        }

        public int getIndex() {
            return index;
        }

        public String getAllele() {
            return allele;
        }
    }
    /**
     * Sorts the chars in the list in ascending order based on the indices
     */
    static String sortByIndices(final String toSort, List<Integer> indices) {
        List<CharWithIndex> charWithIndexList = new ArrayList<>();
        char[] alleles = toSort.toCharArray();
        for (int i = 0; i < indices.size(); i++) {
            charWithIndexList.add(new CharWithIndex(alleles[i], indices.get(i)));
        }
        charWithIndexList.sort(Comparator.comparing(CharWithIndex::getIndex));
        return charWithIndexList.stream().map(CharWithIndex::getAllele).collect(Collectors.joining());
    }

    public int numberOfLabels() {
        return numberOfLabelsPerBase;
    }

}
