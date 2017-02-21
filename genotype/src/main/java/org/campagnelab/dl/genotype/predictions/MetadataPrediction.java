package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * A class to store meta-data (isVariant, isIndel) obtained from the virtual true label.
 * Created by fac2003 on 12/22/16.
 */
public class MetadataPrediction extends Prediction {
    public boolean isVariant;
    public boolean isIndel;
    /**
     * The original goby count index of the reference allele, or -1 if not goby allele/countInfo matched the reference.
     */
    public int referenceGobyIndex;
    /**
     * Permutation, from sorted count indices to original/sbi goby count indices.
     */
    public int[] sorted2OriginalCountIndices;
}
