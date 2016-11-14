package org.campagnelab.dl.varanalysis.learning.mappers;

import org.campagnelab.dl.model.utils.mappers.ConcatFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.NoMaskFeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 5/24/16.
 */
public class ConcatFeatureMapperTest {
    @Test
    public void produceFeature() throws Exception {
        FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> mapper1 = new NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>() {
            @Override
            public int numberOfFeatures() {
                return 4;
            }

            @Override
            public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

            }

            @Override
            public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {

            }

            @Override
            public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
                return featureIndex;
            }
        };
        FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> mapper2 = new NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>() {
            @Override
            public int numberOfFeatures() {
                return 5;
            }

            @Override
            public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

            }

            @Override
            public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {

            }

            @Override
            public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
                return 10 - featureIndex;
            }
        };
        ConcatFeatureMapper concat = new ConcatFeatureMapper(mapper1, mapper2);
        assertEquals(5 + 4, concat.numberOfFeatures());
        float epsilon = 0.01f;
        assertEquals(0f, concat.produceFeature(null, 0), epsilon);
        assertEquals(1f, concat.produceFeature(null, 1), epsilon);
        assertEquals(2f, concat.produceFeature(null, 2), epsilon);
        assertEquals(3f, concat.produceFeature(null, 3), epsilon);

        assertEquals(10f, concat.produceFeature(null, 4), epsilon);
        assertEquals(9f, concat.produceFeature(null, 5), epsilon);
        assertEquals(8f, concat.produceFeature(null, 6), epsilon);
        assertEquals(7f, concat.produceFeature(null, 7), epsilon);
        assertEquals(6f, concat.produceFeature(null, 8), epsilon);


    }

    @Test
    public void mapFeatures() throws Exception {
        FeatureMapper mapper1 = new NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>() {
            @Override
            public int numberOfFeatures() {
                return 4;
            }

            @Override
            public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

            }

            @Override
            public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
                for (int i = 0; i < numberOfFeatures(); i++) {
                    inputs.putScalar(new int[]{i}, produceFeature(record, indexOfRecord));
                }
            }

            @Override
            public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
                return featureIndex;
            }
        };
        FeatureMapper mapper2 = new NoMaskFeatureMapper<BaseInformationRecords.BaseInformationOrBuilder>() {
            @Override
            public int numberOfFeatures() {
                return 5;
            }

            @Override
            public void prepareToNormalize(BaseInformationRecords.BaseInformationOrBuilder record, int indexOfRecord) {

            }

            @Override
            public void mapFeatures(BaseInformationRecords.BaseInformationOrBuilder record, INDArray inputs, int indexOfRecord) {
                for (int i = 0; i < numberOfFeatures(); i++) {
                    inputs.putScalar(new int[]{i}, produceFeature(null, indexOfRecord));
                }
            }

            @Override
            public float produceFeature(BaseInformationRecords.BaseInformationOrBuilder record, int featureIndex) {
                return 10 - featureIndex;
            }
        };
        ConcatFeatureMapper concat = new ConcatFeatureMapper(mapper1, mapper2);
        INDArray inputs = Nd4j.zeros(1, concat.numberOfFeatures());
        concat.prepareToNormalize(null,0);
        concat.mapFeatures(null, inputs, 0);

        assertEquals(5 + 4, concat.numberOfFeatures());
        float epsilon = 0.01f;
        assertEquals(0f, concat.produceFeature(null, 0), epsilon);
        assertEquals(1f, concat.produceFeature(null, 1), epsilon);
        assertEquals(2f, concat.produceFeature(null, 2), epsilon);
        assertEquals(3f, concat.produceFeature(null, 3), epsilon);

        assertEquals(10f, concat.produceFeature(null, 4), epsilon);
        assertEquals(9f, concat.produceFeature(null, 5), epsilon);
        assertEquals(8f, concat.produceFeature(null, 6), epsilon);
        assertEquals(7f, concat.produceFeature(null, 7), epsilon);
        assertEquals(6f, concat.produceFeature(null, 8), epsilon);


    }
}