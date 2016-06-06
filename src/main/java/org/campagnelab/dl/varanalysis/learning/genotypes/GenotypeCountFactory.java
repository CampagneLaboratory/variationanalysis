package org.campagnelab.dl.varanalysis.learning.genotypes;

import org.campagnelab.dl.varanalysis.learning.mappers.GenotypeCount;

/**
 * A factory to encapsulate the creation of spefic types of GenotypeCount implementations.
 * Created by fac2003 on 6/3/16.
 */
public interface GenotypeCountFactory {
    GenotypeCount create();
}
