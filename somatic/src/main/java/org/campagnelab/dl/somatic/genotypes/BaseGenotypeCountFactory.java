package org.campagnelab.dl.somatic.genotypes;

import org.campagnelab.dl.somatic.mappers.GenotypeCount;

/**
 * Created by fac2003 on 6/3/16.
 */
public class BaseGenotypeCountFactory implements GenotypeCountFactory {
    @Override
    public GenotypeCount create() {
        return new GenotypeCount();
    }
}
