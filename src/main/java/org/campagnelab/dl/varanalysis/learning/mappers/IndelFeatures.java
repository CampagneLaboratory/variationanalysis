package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.List;

/**
 * This is the indel feature mapper. maps to binary yes/no isIndel for sorted coutnnts
 *
 * @author Remi Torracinta, rct66
 */

public class IndelFeatures extends AbstractFeatureMapper implements FeatureMapper {

    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:
        return MAX_GENOTYPES * 2;
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
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 : "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES) {
            // germline counts written first:
            final IndelGenotypeCount genotypeCount = (IndelGenotypeCount) getAllCounts(record, false, true).get(featureIndex);
            return genotypeCount.getIsIndel()?1F:0F;
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES;
            final IndelGenotypeCount genotypeCount = (IndelGenotypeCount) getAllCounts(record, true, true).get(featureIndex);
            return genotypeCount.getIsIndel()?1F:0F;
        }
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
}

