package org.campagnelab.dl.model.utils.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.model.utils.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Similar idea to that introduced by Remi in FractionDifferences, but different normalization:
 * here the features are 1/(n+1) where n is the count difference between germline and tumor
 * counts for a given genotype.
 *
 * @author Fabien Campagne
 */

public class FractionDifferences4 extends AbstractFeatureMapper implements FeatureMapper, EfficientFeatureMapper {


    //only implemented for records with 2 samples exactly
    public static final int FRACTION_NORM = 1;
    private int totalCountsGermline;
    private int totalCountsSomatic;


    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:
        // one for each genotype (we ignore strand).
        return AbstractFeatureMapper.MAX_GENOTYPES;
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        totalCountsGermline = 1;
        totalCountsSomatic = 1;
        ObjectArrayList<? extends GenotypeCount> germlineCounts = getAllCounts(record, false, true);
        ObjectArrayList<? extends GenotypeCount> somaticCounts = getAllCounts(record, true, true);
        for (int i = 0; i < germlineCounts.size(); i++) {
            totalCountsGermline += (germlineCounts.get(i).totalCount());
            totalCountsSomatic += (somaticCounts.get(i).totalCount());
        }
        assert(totalCountsGermline>0):"0 total";
        assert(totalCountsSomatic>0):"0 total";
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
        float producedFeat = produceFeatureInternal(record, featureIndex);
        return normalize(producedFeat, FRACTION_NORM);
    }

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, float[] inputs, int offset, int indexOfRecord) {
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            inputs[featureIndex+offset] = produceFeature(record, featureIndex);
        }
    }

    @Override
    public String getFeatureName(int featureIndex) {
        assert featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES : "Only MAX_GENOTYPES features";
        return "fractionDiffCount" + featureIndex;
    }


    private float normalize(float value, int normalizationFactor) {
        assert(value>=0f):"value must be positive " + Float.toString(value);
        assert(value<=1f):"value must be positive " + Float.toString(value);
        return value;
        //return 1f / ((float)(value + 1f));
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert featureIndex >= 0 && featureIndex < AbstractFeatureMapper.MAX_GENOTYPES : "Only MAX_GENOTYPES features";
        ObjectArrayList<? extends GenotypeCount> germlineCounts = getAllCounts(record, false, true);
        ObjectArrayList<? extends GenotypeCount> somaticCounts = getAllCounts(record, true, true);
        // note that we normalize the counts to frequency before substracting:
        int germlineCount = germlineCounts.get(featureIndex).totalCount();
        int somaticCount = somaticCounts.get(featureIndex).totalCount();

        float germFrequency = ((float) germlineCount / (float) totalCountsGermline);
        float somaticFrequency = ((float) somaticCount / (float) totalCountsSomatic);
        // then multiply again by somatic counts to reintroduce size effect:
        return  (float) Math.max(0f,(somaticFrequency - germFrequency));
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {
        return new BaseGenotypeCountFactory() {

            @Override
            public GenotypeCount create() {
                return new GenotypeCount();
            }
        };
    }


    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
    }
}

