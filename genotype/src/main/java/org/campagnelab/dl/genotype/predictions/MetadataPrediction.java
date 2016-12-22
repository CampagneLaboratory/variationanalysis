package org.campagnelab.dl.genotype.predictions;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * A class to store meta-data (isVariant, isIndel) obtained from the virtual true label.
 * Created by fac2003 on 12/22/16.
 */
public class MetadataPrediction extends Prediction {
    public boolean isVariant;
    public boolean isIndel;
}
