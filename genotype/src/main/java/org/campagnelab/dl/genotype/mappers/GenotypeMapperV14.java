package org.campagnelab.dl.genotype.mappers;

/**
 * distinct alleles+inverse+MAX_GENOTYPES=4
 */
public class GenotypeMapperV14 extends GenotypeMapperV11 {
    public GenotypeMapperV14() {
        super();
        this.MAX_GENOTYPES=4;
    }
}
