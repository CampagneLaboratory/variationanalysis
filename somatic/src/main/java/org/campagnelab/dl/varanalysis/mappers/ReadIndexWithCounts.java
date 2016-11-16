package org.campagnelab.dl.varanalysis.mappers;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.List;

/**
 * Created by fac2003 on 6/3/16.
 */
public class ReadIndexWithCounts extends GenotypeCount {
    private IntSet readIndices = new IntOpenHashSet();

    public ReadIndexWithCounts() {
    }

    public void set(List<Integer> readIndicesForwardStrandList, List<Integer> readIndicesReverseStrandList) {
        readIndices.clear();
        readIndices.addAll(readIndicesForwardStrandList);
        readIndices.addAll(readIndicesReverseStrandList);
    }

    public float getDistinctReadIndices() {
        return readIndices.size();
    }
}
