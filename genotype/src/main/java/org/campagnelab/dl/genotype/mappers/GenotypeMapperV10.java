package org.campagnelab.dl.genotype.mappers;

/**
 * Created by rct66 on 12/21/16.
 */
public class GenotypeMapperV10 extends GenotypeMapperV9 {
    public GenotypeMapperV10() {
        super();
        sortCounts = true;
        withDistinctAlleleCounts = false;
        withCombinedLayer = true;
    }

}
