package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.varanalysis.learning.features.Features;
import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapperV2;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 5/25/16.
 */
public class FeatureCalculatorTest {

    /**
     * Check that features generated to DL4J inputs are the same as features generated with produceFeature:
     */
    @Test
    public void makeCheckFeatures() {

        FeatureMapperV2 v2 = new FeatureMapperV2();
        BaseInformationRecords.BaseInformation record = prepareRecord();
        v2.prepareToNormalize(record, 0);
        Features features = new Features(v2.numberOfFeatures());
        for (int i = 0; i < v2.numberOfFeatures(); i++) {
            features.setFeatureValue(v2.produceFeature(record, i), i);
        }
        INDArray inputs = Nd4j.zeros(1, v2.numberOfFeatures());
        v2.mapFeatures(record, inputs, 0);
        Features featuresFromArray=new Features(inputs,1);
        assertEquals(features, featuresFromArray);
        System.out.println(features);
    }

    private BaseInformationRecords.BaseInformation prepareRecord() {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setPosition(1);
        builder.setReferenceIndex(0);
        builder.setReferenceBase("A");


        //germline counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilder = BaseInformationRecords.SampleInfo.newBuilder();
        BaseInformationRecords.CountInfo.Builder builderInfo = BaseInformationRecords.CountInfo.newBuilder();
        builderInfo.setFromSequence("A");
        builderInfo.setToSequence("C");
        builderInfo.setMatchesReference(true);
        builderInfo.setGenotypeCountForwardStrand(10);
        builderInfo.setGenotypeCountReverseStrand(20);
        sampleBuilder.addCounts(builderInfo.build());
        builder.addSamples(sampleBuilder.build());

        //somatic counts
        BaseInformationRecords.SampleInfo.Builder sampleBuilderS = BaseInformationRecords.SampleInfo.newBuilder();
        sampleBuilderS.setIsTumor(true);
        BaseInformationRecords.CountInfo.Builder builderInfoS = BaseInformationRecords.CountInfo.newBuilder();
        builderInfoS.setFromSequence("A");
        builderInfoS.setToSequence("T");
        builderInfoS.setMatchesReference(true);
        builderInfoS.setGenotypeCountForwardStrand(50);
        builderInfoS.setGenotypeCountReverseStrand(50);
        sampleBuilderS.addCounts(builderInfoS.build());
        builder.addSamples(sampleBuilderS.build());
        return builder.build();
    }
}
