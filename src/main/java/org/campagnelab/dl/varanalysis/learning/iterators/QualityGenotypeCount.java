package org.campagnelab.dl.varanalysis.learning.iterators;

/**
 * Created by rct66 on 6/1/16.
 */
public class QualityGenotypeCount extends GenotypeCount {

    int forwardCount;
    int reverseCount;
    String toSequence;
    float qualityScoreForward;
    float qualityScoreReverse;

    public QualityGenotypeCount(int forwardCount, int reverseCount, String toSequence, float qualityScoreForward, float qualityScoreReverse) {
        super(forwardCount, reverseCount, toSequence);
        this.qualityScoreForward = qualityScoreForward;
        this.qualityScoreReverse = qualityScoreReverse;
    }


    public int totalCount() {
        return forwardCount + reverseCount;
    }

    public float getQualityScoreForward(){
        return qualityScoreForward;
    }

    public float getQualityScoreReverse(){
        return qualityScoreReverse;
    }

    @Override
    public int compareTo(GenotypeCount o) {
        return o.totalCount() - totalCount();
    }

    @Override
    public String toString() {
        return String.format("totalCount=%d %d on + / %d on - %s, quality on + / %e on - %e",totalCount(), forwardCount,reverseCount,toSequence,qualityScoreForward,qualityScoreReverse);
    }

}
