package org.campagnelab.dl.framework.mappers.processing;

import org.campagnelab.dl.framework.mappers.TwoDimensionalConcatFeatureMapper;
import org.campagnelab.dl.framework.mappers.RNNFeatureMapper;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 12/13/16.
 */
public class TwoDimensionalRemoveMaskFeatureMapperTest {
    @Test
    public void testFeatureFilterConcatRNNFeatureMapperOHBMDelegates() {
        String sequence1 = "ATC";
        String sequence2 = "GN";
        String sequence3 = "J";
        String expectedFeatures =
                "[[[1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]]]";
        String expectedMask =
                "[[1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],\n" +
                " [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.00],\n" +
                " [1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]]";
        RNNFeatureMapper<String> rnnFeatureMapper1 = new RNNFeatureMapper<>(3, Function.identity(), String::length);
        RNNFeatureMapper<String> rnnFeatureMapper2 = new RNNFeatureMapper<>(3, Function.identity(), String::length);
        RNNFeatureMapper<String> rnnFeatureMapper3 = new RNNFeatureMapper<>(3, Function.identity(), String::length);
        TwoDimensionalConcatFeatureMapper<String> twoDConcatFMapper = new TwoDimensionalConcatFeatureMapper<>(rnnFeatureMapper1, rnnFeatureMapper2, rnnFeatureMapper3);
        TwoDimensionalRemoveMaskFeatureMapper<String> twoDFilterFMapper = new TwoDimensionalRemoveMaskFeatureMapper<>(twoDConcatFMapper);
        twoDFilterFMapper.prepareToNormalize(sequence1, 0);
        twoDFilterFMapper.prepareToNormalize(sequence2, 1);
        twoDFilterFMapper.prepareToNormalize(sequence3, 2);
        INDArray inputs = Nd4j.zeros(3, 6, 9);
        twoDFilterFMapper.mapFeatures(sequence1, inputs, 0);
        twoDFilterFMapper.mapFeatures(sequence2, inputs, 1);
        twoDFilterFMapper.mapFeatures(sequence3, inputs, 2);
        assertEquals(inputs.toString(), expectedFeatures);
        INDArray mask = Nd4j.zeros(3, 9);
        twoDFilterFMapper.maskFeatures(sequence1, mask, 0);
        twoDFilterFMapper.maskFeatures(sequence2, mask, 1);
        twoDFilterFMapper.maskFeatures(sequence3, mask, 2);
        assertEquals(mask.toString(), expectedMask);
        INDArray sequence1Features = inputs.getRow(0);
        INDArray sequence2Mask = mask.getRow(1);
        for (int i = 0; i < 9; i++) {
            assertEquals(sequence2Mask.getInt(i) == 1, twoDFilterFMapper.isMasked(sequence2, i * 6));
            for (int j = 0; j < 6; j++) {
                int featureIndex = i * 6 + j;
                assertEquals(sequence1Features.getFloat(j, i), twoDFilterFMapper.produceFeature(sequence1, featureIndex), 1e-9);
            }
        }
    }
}
