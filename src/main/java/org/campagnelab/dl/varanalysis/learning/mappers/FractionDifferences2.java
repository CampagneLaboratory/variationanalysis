package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Similar idea to that introduced by Remi in FractionDifferences, but different normalization:
 * here the features are 1/(n+1) where n is the count difference between germline and tumor
 * counts for a given genotype.
 *
 * @author Fabien Campagne
 */

public class FractionDifferences2 extends AbstractFeatureMapper implements FeatureMapper {


    //only implemented for records with 2 samples exactly
    public static final int FRACTION_NORM = 1;
    private int totalCountsGermline;
    private int totalCountsSomatic;


    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:
        // one for each genotype (we ignore strand).
        return MAX_GENOTYPES;
    }

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        totalCountsGermline = 0;
        totalCountsSomatic = 0;
        for (int i = 0; i < MAX_GENOTYPES; i++) {
            BaseInformationRecords.CountInfo germline = record.getSamples(0).getCounts(i);
            BaseInformationRecords.CountInfo somatic = record.getSamples(1).getCounts(i);
            totalCountsGermline += (germline.getGenotypeCountForwardStrand() + germline.getGenotypeCountReverseStrand());
            totalCountsSomatic += (somatic.getGenotypeCountForwardStrand() + somatic.getGenotypeCountReverseStrand());
        }
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
        return normalize(produceFeatureInternal(record, featureIndex), FRACTION_NORM);
    }


    private float normalize(float value, int normalizationFactor) {
        return 1f / ((float)(value + 1f));
    }


    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {

        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES : "Only MAX_GENOTYPES features";
        ObjectArrayList<? extends GenotypeCount> germlineCounts = getAllCounts(record, false, false);
        ObjectArrayList<? extends GenotypeCount> somaticCounts = getAllCounts(record, true, false);
        // note that we normalize the counts to frequency before substracting:
        int germlineCount = germlineCounts.get(featureIndex).forwardCount + germlineCounts.get(featureIndex).reverseCount;
        int somaticCount = somaticCounts.get(featureIndex).forwardCount + somaticCounts.get(featureIndex).reverseCount;

        float germFrequency = ((float) germlineCount / totalCountsGermline);
        float somaticFrequency = ((float) somaticCount) / totalCountsSomatic;
        // then multiply again by somatic counts to reintroduce size effect:
        return  (float) (Math.max(0f,(somaticFrequency - germFrequency)) * somaticCount);
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

