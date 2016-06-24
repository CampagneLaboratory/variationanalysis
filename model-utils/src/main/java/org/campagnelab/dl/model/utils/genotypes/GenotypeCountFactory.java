package org.campagnelab.dl.model.utils.genotypes;

import org.campagnelab.dl.model.utils.mappers.GenotypeCount;

/**
 * A factory to encapsulate the creation of spefic types of GenotypeCount implementations.
 * Created by fac2003 on 6/3/16.
 */
public interface GenotypeCountFactory {
    GenotypeCount create();
}
