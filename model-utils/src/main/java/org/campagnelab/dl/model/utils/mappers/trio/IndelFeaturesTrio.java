package org.campagnelab.dl.model.utils.mappers.trio;

import org.campagnelab.dl.model.utils.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.model.utils.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.model.utils.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.GenotypeCount;
import org.campagnelab.dl.model.utils.mappers.IndelGenotypeCount;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * This is the indel feature mapper. maps to binary yes/no isIndel for sorted coutnnts
 *
 * @author Remi Torracinta, rct66
 */

public class IndelFeaturesTrio extends AbstractFeatureMapperTrio<BaseInformationRecords.BaseInformationOrBuilder>
        implements FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> {

    public int numberOfFeatures() {
        // now including features for both parent samples
        return AbstractFeatureMapper.MAX_GENOTYPES;
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //shouldn't need to do anything
    }


    int[] indices = new int[]{0, 0};

    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return produceFeatureInternal(record, featureIndex);
    }


    private float normalize(float value, int normalizationFactor) {
        if (normalizationFactor == 0) {
            return 0;
        }
        float normalized = value / normalizationFactor;
        assert normalized >= 0 && normalized <= 1 : "value must be normalized: " + normalized;
    return normalized;
}


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES : "Only MAX_GENOTYPES*2 features";
        // father counts written first
        final IndelGenotypeCount genotypeCount = (IndelGenotypeCount) getAllCounts(record, 0, true).get(featureIndex);
        return genotypeCount.getIsIndel() ? 1F : 0F;

    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {
        return new BaseGenotypeCountFactory() {

            @Override
            public GenotypeCount create() {
                return new IndelGenotypeCount();
            }
        };
    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        IndelGenotypeCount myCount = (IndelGenotypeCount) count;
        myCount.set(sampleCounts.getIsIndel());
    }

    @Override
    public String getFeatureName(int featureIndex) {
        assert featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES : "Only MAX_GENOTYPES features";
        return "count" + featureIndex + "isIndel";
    }
}

