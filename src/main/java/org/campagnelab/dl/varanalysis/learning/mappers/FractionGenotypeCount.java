package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * Created by rct66 on 6/1/16.
 */
public class FractionGenotypeCount extends GenotypeCount {

    /**
     * Fraction genotypecount subclass for the unused fractiondifference class. Not complete.
     */

    private float fractionDifference;

    public FractionGenotypeCount() {
    }

    public float getFractionDifference() {
        return fractionDifference;
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s, quality on + / %e on - %e", totalCount(), forwardCount, reverseCount, toSequence, fractionDifference);
    }

    public void set(float fractionDifference) {
        this.fractionDifference = fractionDifference;
    }
}
