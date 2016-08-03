package org.campagnelab.dl.model.utils.mappers;

import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.model.utils.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.model.utils.mappers.AbstractFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.GenotypeCount;
import org.campagnelab.dl.model.utils.mappers.ReadIndexWithCounts;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Estimates the number of distinct read indices per base position.
 * Created by fac2003 on 6/3/16.
 *
 * @author Fabien Campagne
 */
public class ReadIndexFeaturesFix extends AbstractFeatureMapper implements FeatureMapper {

    public static final int NUM_GENOTYPES = 5;
    public static final int NUM_SAMPLES = 2;

    @Override
    public int numberOfFeatures() {
        // only one feature per genotype: the number of distinct read indices in each sample. 5 genotypes max.
        return NUM_SAMPLES * NUM_GENOTYPES;
    }

    public static final int MAX_GENOTYPES = 5;

    int sumReadIndex;

    /**
     *  normalization is ndone by taking the inverse of the number of distinct read indices.
     * @param record        The record to convert to features & labels.
     * @param indexOfRecord Index of the record in the destination dataset.
     */
    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

//        indices[0] = indexOfRecord;
//        sumReadIndex = 0;
//        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
//            sumReadIndex += produceFeatureInternal(record, featureIndex);
//        }
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        prepareToNormalize(record, indexOfRecord);
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }

    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        return normalize(produceFeatureInternal(record, featureIndex), sumReadIndex);
    }

    @Override
    public String getFeatureName(int featureIndex) {
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2: "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES) {
            return "numReadIdxsGermlineCount"+(featureIndex/2);
        } else {
            return "numReadIdxsSomaticCount"+(featureIndex/2);
        }
    }


    private float normalize(float value, int normalizationFactor) {
        // produce a value between zero and 1, increasingly smaller as the number of distinct read indices increases.
        // plus one is to prevent division by zero for genotypes without any read.
        return 1f / (value+1);
        /*if (normalizationFactor == 0) {
            return 0;
        }
        float normalized = value / normalizationFactor;
        assert normalized >= 0 && normalized <= 1 : "value must be normalized: " + normalized;
        return normalized;*/
    }

    //TODO: at the moment, every readindex count from the germline sample is repeated twice, and those of somatic are missing.
    public float produceFeatureInternal(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2: "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES) {
            // germline counts written first:
            final ReadIndexWithCounts genotypeCount = (ReadIndexWithCounts) getAllCounts(record, false).get(featureIndex / 2);
            return genotypeCount.getDistinctReadIndices();
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES;
            final ReadIndexWithCounts genotypeCount = (ReadIndexWithCounts) getAllCounts(record, true).get(featureIndex / 2);
            return genotypeCount.getDistinctReadIndices();
        }
    }

    public boolean oneSampleHasTumor(java.util.List<BaseInformationRecords.SampleInfo> samples) {
        for (BaseInformationRecords.SampleInfo sample : samples) {
            if (sample.getIsTumor()) return true;
        }
        return false;

    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        ReadIndexWithCounts myCounts = (ReadIndexWithCounts) count;
        myCounts.set(ProtoPredictor.expandFreq(sampleCounts.getReadIndicesForwardStrandList()),
                ProtoPredictor.expandFreq(sampleCounts.getReadIndicesReverseStrandList()));
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory() {
            @Override
            public GenotypeCount create() {
                return new ReadIndexWithCounts();
            }
        };
    }

}
