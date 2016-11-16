package org.campagnelab.dl.somatic.mappers;

/**
 * Another feature mapper. Based of V9, but with FractionDifferences4 fix when counts==0;
 */
public class FeatureMapperV9_2 extends ConcatFeatureMapper {
    public FeatureMapperV9_2() {
        super(new SimpleFeatureCalculator(true),
                new SimpleFeatureCalculator(false),
                new SortedGenotypeAgreementMapper(),
                //new MagnitudeFeatures(),
                // new QualityFeatures(),
                new ReadIndexFeatures(),
                new FractionDifferences4()

        );
    }
}
