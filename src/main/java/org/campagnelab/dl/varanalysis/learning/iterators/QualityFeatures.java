package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collections;

/**
 * This is the quality score feature mapper. maps phred averages from parquet
 *
 * @author Remi Torracinta, rct66
 */

public class QualityFeatures implements FeatureCalculator {


    public static final int MAX_GENOTYPES = 5;

    @Override
    public int numberOfFeatures() {
        // we need features for the normal sample and for the tumor sample:

        return MAX_GENOTYPES * 2*2;
    }

    int sumCounts;

    public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {
        //shouldn't need to do anything
    }

    @Override
    public int numberOfLabels() {
        return 20;
    }

    int[] indices = new int[]{0, 0};

    @Override
    public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
        indices[0] = indexOfRecord;
        for (int featureIndex = 0; featureIndex < numberOfFeatures(); featureIndex++) {
            indices[1] = featureIndex;
            inputs.putScalar(indices, produceFeature(record, featureIndex));
        }
    }


    @Override
    public void mapLabels(BaseInformationRecords.BaseInformationOrBuilder record, INDArray labels, int indexOfRecord) {
        indices[0] = indexOfRecord;

        for (int labelIndex = 0; labelIndex < numberOfLabels(); labelIndex++) {
            indices[1] = labelIndex;
            labels.putScalar(indices, produceLabel(record, labelIndex));
        }
    }


    public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
        assert featureIndex >= 0 && featureIndex < MAX_GENOTYPES * 2*2 : "Only MAX_GENOTYPES*2*2 features";
        if (featureIndex < MAX_GENOTYPES * 2) {
            // germline counts written first:
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return record.getSamples(0).getCounts(featureIndex/2).getQualityScoreForwardStrand();
            } else {
                return record.getSamples(0).getCounts(featureIndex/2).getQualityScoreReverseStrand();
            }
        } else {
            // tumor counts written next:
            featureIndex -= MAX_GENOTYPES * 2;
            if ((featureIndex % 2) == 1) {
                // odd featureIndices are forward strand:
                return record.getSamples(1).getCounts(featureIndex/2).getQualityScoreForwardStrand();
            } else {
                return record.getSamples(1).getCounts(featureIndex/2).getQualityScoreReverseStrand();
            }
        }
    }


    }
}
