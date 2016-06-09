package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.learning.genotypes.BaseGenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.genotypes.GenotypeCountFactory;
import org.campagnelab.dl.varanalysis.learning.iterators.AbstractFeatureMapper;
import org.campagnelab.dl.varanalysis.learning.iterators.FeatureCalculator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collections;

/**
 * This is a simple feature mapper. It is designed using information currently available in the Parquet file.
 * Each position has the following information:
 * <pre> 69033640	11	false
 * position=14521   referenceIndex=0       isMutated=false
 * sample 0 counts=[10, 0, 0, 63, 0, 4, 0, 1, 62, 0]
 * sample 1 counts=[2, 0, 0, 45, 0, 3, 0, 0, 68, 0]
 *
 * position=14521   referenceIndex=0       isMutated=true
 * sample 0 counts=[10, 0, 0, 63, 0, 4, 0, 1, 62, 0]
 * sample 1 counts=[2, 11, 0, 34, 0, 3, 12, 0, 56, 0]
 * </pre>
 * <p>
 * Using these data, we can map to features as follows:
 * <ul>
 * <li>Make the isMutated boolean the only label. This is not ideal: the net may tell us if the base is mutated,
 * but we will not know what the mutation is..</li>
 * <li>Concatenate the count integers and use these as features. The only way for the net to learn from these data is to count
 * the number of counts elements that has "enough" reads to call a genotype. If the genotype calls are more than two, then
 * the site is likely mutated, because most sites will be heterozygous at most. With these features, I expect the
 * net to have more trouble predicting mutations at homozygous sites, than at heterozygous sites. We'll see. </ul>
 * Created by fac2003 on 5/21/16.
 *
 * @author Fabien Campagne
 */

public class SimpleFeatureCalculator extends AbstractFeatureMapper implements FeatureCalculator {


    public SimpleFeatureCalculator(boolean sort){
        this.sort = sort;
    }

    public SimpleFeatureCalculator(){
        this.sort = true;
    }

    @Override
    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:

        return MAX_GENOTYPES * 2 * 2;
    }
    boolean sort; //by default we sort
    int sumCounts;

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        indices[0] = indexOfRecord;
        sumCounts = 0;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            sumCounts += produceFeatureInternal(record, featureIndex);
        }
    }

    @Override
    public int numberOfLabels() {
        return 2;
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
        return normalize(produceFeatureInternal(record, featureIndex), sumCounts);
    }

    @Override
    public void mapLabels(BaseInformationRecords.BaseInformationOrBuilder record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
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
        assert (featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2 * 2) : "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, false, sort).get(featureIndex / 2).forwardCount;
            } else {
                return getAllCounts(record, false, sort).get(featureIndex / 2).reverseCount;
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return getAllCounts(record, true, sort).get(featureIndex / 2).forwardCount;
            } else {
                return getAllCounts(record, true, sort).get(featureIndex / 2).reverseCount;
            }
        }
    }


    @Override
    public float produceLabel(BaseInformationRecords.BaseInformationOrBuilder record, int labelIndex) {
        assert labelIndex == 0 || labelIndex == 1 : "only one label.";

        //return record.getMutated() ? 1.0f : 0.0f;
        if (labelIndex == 0) return record.getMutated() ? 1 : 0;
        else {
            return !record.getMutated() ? 1 : 0;
        }

    }

    @Override
    protected void initializeCount(BaseInformationRecords.CountInfo sampleCounts, GenotypeCount count) {
        // nothing to do, already done in the base class.
    }

    @Override
    protected GenotypeCountFactory getGenotypeCountFactory() {

        return new BaseGenotypeCountFactory();

    }
}
