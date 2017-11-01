package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

import java.util.BitSet;

/**
 * Stores meta data for segments: whether each base matches ref, a SNP or an indel.
 * Created by fac2003 on 11/1/17.
 */
public class SegmentMetaData extends Prediction {
    public SegmentMetaData() {
        isSnp = new BitSet();
        isIndel = new BitSet();
    }

    public BitSet isSnp;
    public BitSet isIndel;
}
