package org.campagnelab.dl.varanalysis.mappers;

/**
 * Created by fac2003 on 5/25/16.
 */
public class GenotypeCount implements Comparable<GenotypeCount> {
    public int forwardCount;
    public int reverseCount;
    int compareCount;
    String toSequence;
    int genotypeIndex;

    public GenotypeCount() {
    }

    public GenotypeCount(int forwardCount, int reverseCount, String toSequence, int genotypeIndex,int compareCount) {
        set(forwardCount, reverseCount, toSequence,genotypeIndex, compareCount);
    }

    public void set(int forwardCount, int reverseCount, String toSequence,int genotypeIndex,int compareCount) {
        this.forwardCount = forwardCount;
        this.reverseCount = reverseCount;
        this.toSequence = toSequence;
        this.genotypeIndex=genotypeIndex;
        this.compareCount = compareCount;
    }

    public int totalCount() {
        return forwardCount + reverseCount;
    }

    public int getCompareCount() {
        return compareCount;
    }

    @Override
    public int compareTo(GenotypeCount o) {
        return o.getCompareCount() - getCompareCount();
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s", totalCount(), forwardCount, reverseCount, toSequence);
    }
}
