package org.campagnelab.dl.varanalysis.learning.mappers;

/**
 * Created by rct66 on 6/1/16.
 */
public class QualityGenotypeCount extends GenotypeCount {


    private float qualityScoreForward;
    private float qualityScoreReverse;

    public QualityGenotypeCount() {
    }


    public int totalCount() {
        return forwardCount + reverseCount;
    }

    public float getQualityScoreForward() {
        return qualityScoreForward;
    }

    public float getQualityScoreReverse() {
        return qualityScoreReverse;
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s, quality on + / %e on - %e", totalCount(), forwardCount, reverseCount, toSequence, qualityScoreForward, qualityScoreReverse);
    }

    public void set(float averageQualityForward, float averageQualityReverse) {
        this.qualityScoreForward = averageQualityForward;
        this.qualityScoreReverse = averageQualityReverse;
    }
}
