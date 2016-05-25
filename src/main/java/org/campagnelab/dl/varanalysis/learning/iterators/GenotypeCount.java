package org.campagnelab.dl.varanalysis.learning.iterators;

/**
 * Created by fac2003 on 5/25/16.
 */
public class GenotypeCount implements Comparable<GenotypeCount> {
    int forwardCount;
    int reverseCount;
    String toSequence;

    public GenotypeCount(int forwardCount, int reverseCount, String toSequence) {
        this.forwardCount = forwardCount;
        this.reverseCount = reverseCount;
        this.toSequence = toSequence;
    }

    public int totalCount() {
        return forwardCount + reverseCount;
    }

    @Override
    public int compareTo(GenotypeCount o) {
        return o.totalCount() - totalCount();
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s",totalCount(), forwardCount,reverseCount,toSequence);
    }
}
