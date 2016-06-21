package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * Created by rct66 on 6/1/16.
 */
public class IndelGenotypeCount extends GenotypeCount {


    private boolean isIndel;

    public IndelGenotypeCount() {
    }

    public boolean getIsIndel() {
        return isIndel;
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s, isIndel %b", totalCount(), forwardCount, reverseCount, toSequence, isIndel);
    }

    public void set(boolean isIndel) {
        this.isIndel = isIndel;
    }
}
