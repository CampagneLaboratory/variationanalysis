package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.List;

/**
 * Created by fac2003 on 6/3/16.
 */
public class ReadIndexWithCounts extends GenotypeCount {
    IntSet readIndices = new IntArraySet();

    public ReadIndexWithCounts(int forwardCount, int reverseCount, String toSequence, List<Integer>
            readIndicesForwardStrandList, List<Integer> readIndicesReverseStrandList) {
        super(forwardCount, reverseCount, toSequence);
        readIndices.addAll(readIndicesForwardStrandList);
        readIndices.addAll(readIndicesReverseStrandList);
    }

    public float getDistinctReadIndices() {
        return readIndices.size();
    }
}
