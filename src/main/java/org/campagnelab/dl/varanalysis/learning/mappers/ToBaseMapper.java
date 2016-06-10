package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
/**
 * Feaature mapper for option 2:
 * Option 2.  Add boolean to indicate when the bases match: 1 or don’t (0).
 * <p>
 * e.g.,
 * [100,100,0,0]
 * [100,100,0,0] 1, 0, 0, 1 .
 * Problem with ties on the counts, will produce 0 when 1 would be valid.
 * <p>
 * Created by fac2003 on 6/10/16.
 */
public class ToBaseMapper extends AbstractFeatureMapper implements FeatureMapper {

    private static final int MAX_GENOTYPES = 5;
    private int[] indices = new int[1];


    @Override
    public int numberOfFeatures() {
        return MAX_GENOTYPES;
    }

    @Override
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {

        indices[0] = indexOfRecord;
        prepareToNormalize(record, indexOfRecord);
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }


    }

    @Override
    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
// compare genotype indices after sorting:
        int genotypeIndexSample0 = getAllCounts(record.getSamples(0), getGenotypeCountFactory(), true).get(featureIndex).genotypeIndex;
        int genotypeIndexSample1 = getAllCounts(record.getSamples(0), getGenotypeCountFactory(), true).get(featureIndex).genotypeIndex;
        return genotypeIndexSample0 == genotypeIndexSample1 ? 1f : .0f;
    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        // do nothing. We already have the genotype index in the count.
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {
        return new BaseGenotypeCountFactory();
    }
}
