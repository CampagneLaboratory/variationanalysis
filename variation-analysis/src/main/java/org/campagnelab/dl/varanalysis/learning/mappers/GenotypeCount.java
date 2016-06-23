package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * Created by fac2003 on 5/25/16.
 */
public class GenotypeCount implements Comparable<GenotypeCount> {
    int forwardCount;
    int reverseCount;
    int germCount;
    String toSequence;
    int genotypeIndex;

    public GenotypeCount() {
    }

    public GenotypeCount(int forwardCount, int reverseCount, String toSequence, int genotypeIndex,int germCount) {
        set(forwardCount, reverseCount, toSequence,genotypeIndex, germCount);
    }

    public void set(int forwardCount, int reverseCount, String toSequence,int genotypeIndex,int germCount) {
        this.forwardCount = forwardCount;
        this.reverseCount = reverseCount;
        this.toSequence = toSequence;
        this.genotypeIndex=genotypeIndex;
        this.germCount = germCount;
    }

    public int totalCount() {
        return forwardCount + reverseCount;
    }

    public int totalCountGerm() {
        return germCount;
    }

    @Override
    public int compareTo(GenotypeCount o) {
        return o.totalCountGerm() - totalCountGerm();
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s", totalCount(), forwardCount, reverseCount, toSequence);
    }
}
