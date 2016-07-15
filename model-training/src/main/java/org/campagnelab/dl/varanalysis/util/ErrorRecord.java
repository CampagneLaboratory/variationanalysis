package org.campagnelab.dl.varanalysis.util;

import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Comparator;

/**
 * Created by fac2003 on 7/15/16.
 */
public class ErrorRecord {
    /* A measure of how wrong this prediction was. Larger is more wronger.
         */
    float wrongness;
    INDArray features;

    INDArray label;

    public ErrorRecord(float wrongness, INDArray features, INDArray label) {
        this.wrongness = wrongness;
        this.features = features;
        this.label = label;
    }

    public static Comparator<ErrorRecord> INCREASING_SCORE_COMPARATOR = new Comparator<ErrorRecord>() {
        public int compare(ErrorRecord a, ErrorRecord b) {
            return Float.compare(a.wrongness, b.wrongness);
        }
    };

    @Override
    public String toString() {
        return String.format("%f features=%s labels=%s",wrongness,features,label);
    }
}
