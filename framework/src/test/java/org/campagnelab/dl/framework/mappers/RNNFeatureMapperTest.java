package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;


/**
 * Created by joc2080 on 11/28/16.
 */
public class RNNFeatureMapperTest {
    @Test
    public void mapRNNFeaturesOHBFMDelegates() {
        String sequence1 = "ATCGNJ";
        String sequence2 = "ATCG";
        String sequence3 = "AT";
        String expectedFeatures =
                "[[[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 1.00]],\n" +
                "\n" +
                " [[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]]]";
        String expectedMask =
                "[[1.00, 1.00, 1.00, 1.00, 1.00, 1.00],\n" +
                " [1.00, 1.00, 1.00, 1.00, 0.00, 0.00],\n" +
                " [1.00, 1.00, 0.00, 0.00, 0.00, 0.00]]";
        RNNFeatureMapper<String> rnnFeatureMapper = new RNNFeatureMapper<>(6, Function.identity(), String::length);

        INDArray inputs = Nd4j.zeros(3, 6, 6);
        INDArray mask = Nd4j.zeros(3, 6);

        rnnFeatureMapper.prepareToNormalize(sequence1, 0);
        rnnFeatureMapper.mapFeatures(sequence1, inputs, 0);
        rnnFeatureMapper.maskFeatures(sequence1, mask, 0);

        INDArray sequence1Features = inputs.getRow(0);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                int featureIndex = i * 6 + j;
                assertEquals(sequence1Features.getFloat(j, i), rnnFeatureMapper.produceFeature(sequence1, featureIndex), 1e-9);
            }
        }

        rnnFeatureMapper.prepareToNormalize(sequence2, 1);
        rnnFeatureMapper.mapFeatures(sequence2, inputs, 1);
        rnnFeatureMapper.maskFeatures(sequence2, mask, 1);

        INDArray sequence2Mask = mask.getRow(1);
        for (int i = 0; i < 6; i++) {
            assertEquals(sequence2Mask.getInt(i) == 1, rnnFeatureMapper.isMasked(sequence2, i * 6));
        }

        rnnFeatureMapper.prepareToNormalize(sequence3, 2);
        rnnFeatureMapper.mapFeatures(sequence3, inputs, 2);
        rnnFeatureMapper.maskFeatures(sequence3, mask, 2);

        assertEquals(inputs.toString(), expectedFeatures);
        assertEquals(mask.toString(), expectedMask);
    }
}
