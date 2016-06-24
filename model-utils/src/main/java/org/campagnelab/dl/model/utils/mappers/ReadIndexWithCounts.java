package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.List;

/**
 * Created by fac2003 on 6/3/16.
 */
public class ReadIndexWithCounts extends GenotypeCount {
    private IntSet readIndices = new IntArraySet();

    public ReadIndexWithCounts() {
    }

    public void set(List<Integer>
            readIndicesForwardStrandList, List<Integer> readIndicesReverseStrandList) {

        readIndices.addAll(readIndicesForwardStrandList);
        readIndices.addAll(readIndicesReverseStrandList);
    }

    public float getDistinctReadIndices() {
        return readIndices.size();
    }
}
