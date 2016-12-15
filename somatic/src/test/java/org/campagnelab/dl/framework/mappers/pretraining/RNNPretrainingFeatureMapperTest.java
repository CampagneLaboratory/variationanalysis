package org.campagnelab.dl.framework.mappers.pretraining;

import org.campagnelab.dl.framework.mappers.RNNFeatureMapper;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

/**
 * Created by joshuacohen on 12/15/16.
 */
public class RNNPretrainingFeatureMapperTest {
    @Test
    public void testRNNPretrainingMapperNoEOS() {
        String sequence1 = "ATC";
        String sequence2 = "GN";
        String sequence3 = "J";
        String expectedFeatures =
                "[[[1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00]]]";
        String expectedMask =
                "[[1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],\n" +
                " [1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00],\n" +
                " [1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00]]";
        RNNFeatureMapper<String> rnnFeatureMapper = new RNNFeatureMapper<>(3, Function.identity(),
                String::length);
        INDArray inputs1 = Nd4j.zeros(3, 6, 3);
        rnnFeatureMapper.mapFeatures(sequence1, inputs1, 0);
        rnnFeatureMapper.mapFeatures(sequence2, inputs1, 1);
        rnnFeatureMapper.mapFeatures(sequence3, inputs1, 2);
        RNNPretrainingFeatureMapper<String> rnnPretrainingFeatureMapper = new RNNPretrainingFeatureMapper<>(
                rnnFeatureMapper, null, String::length);
        rnnPretrainingFeatureMapper.prepareToNormalize(sequence1, 0);
        rnnPretrainingFeatureMapper.prepareToNormalize(sequence2, 1);
        rnnPretrainingFeatureMapper.prepareToNormalize(sequence3, 2);
        INDArray inputs = Nd4j.zeros(3, 7, 7);
        rnnPretrainingFeatureMapper.mapFeatures(sequence1, inputs, 0);
        rnnPretrainingFeatureMapper.mapFeatures(sequence2, inputs, 1);
        rnnPretrainingFeatureMapper.mapFeatures(sequence3, inputs, 2);
        INDArray mask = Nd4j.zeros(3, 7);
        rnnPretrainingFeatureMapper.maskFeatures(sequence1, mask, 0);
        rnnPretrainingFeatureMapper.maskFeatures(sequence2, mask, 1);
        rnnPretrainingFeatureMapper.maskFeatures(sequence3, mask, 2);
    }
}
