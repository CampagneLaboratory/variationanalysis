package org.campagnelab.dl.genotype.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.Prediction;

/**
 * Represent the number of distint alleles called at a site.
 *  Created by fac2003 on 12/20/16.
 */
public class NumDistinctAlleles extends Prediction {
    int trueValue;
    int predictedValue;
}
