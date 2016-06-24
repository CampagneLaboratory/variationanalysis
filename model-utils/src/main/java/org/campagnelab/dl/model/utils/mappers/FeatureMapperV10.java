package org.campagnelab.dl.model.utils.mappers;

/**
 * Another feature mapper.
 */
public class FeatureMapperV10 extends ConcatFeatureMapper {
    public FeatureMapperV10() {
        super(new SimpleFeatureCalculator(true),
                new SimpleFeatureCalculator(false),
                new SortedGenotypeAgreementMapper(),
                //new MagnitudeFeatures(),
                new QualityFeatures(),
                new ReadIndexFeatures(),
                new FractionDifferences2()

        );
    }
}
