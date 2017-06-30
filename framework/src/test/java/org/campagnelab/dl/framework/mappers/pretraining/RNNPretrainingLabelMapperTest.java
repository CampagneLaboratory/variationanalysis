package org.campagnelab.dl.framework.mappers.pretraining;

import org.campagnelab.dl.framework.mappers.RNNFeatureMapper;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

import static org.junit.Assert.assertEquals;


/**
 * Created by joshuacohen on 12/15/16.
 */
public class RNNPretrainingLabelMapperTest {
    @Test
    public void testRNNPretrainingLabelMapperNoEOS() {
        String sequence1 = "ATC";
        String sequence2 = "GN";
        String sequence3 = "J";
        String expectedLabels =
                "[[[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00]],\n" +
                "\n" +
                " [[0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00]]]";
        String expectedMask =
                "[[0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00],\n" +
                " [0.00, 0.00, 1.00, 1.00, 1.00, 0.00, 0.00],\n" +
                " [0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00]]";
        RNNFeatureMapper<String> rnnFeatureMapper = new RNNFeatureMapper<>(3, Function.identity(),
                String::length);
        RNNPretrainingLabelMapper<String> rnnPretrainingLabelMapper = new RNNPretrainingLabelMapper<>(
                rnnFeatureMapper, null, String::length);
        INDArray inputs = Nd4j.zeros(3, 8, 7);
        INDArray mask = Nd4j.zeros(3, 7);

        rnnPretrainingLabelMapper.prepareToNormalize(sequence1, 0);
        rnnPretrainingLabelMapper.mapLabels(sequence1, inputs, 0);
        rnnPretrainingLabelMapper.maskLabels(sequence1, mask, 0);

        rnnPretrainingLabelMapper.prepareToNormalize(sequence2, 1);
        rnnPretrainingLabelMapper.mapLabels(sequence2, inputs, 1);
        rnnPretrainingLabelMapper.maskLabels(sequence2, mask, 1);

        rnnPretrainingLabelMapper.prepareToNormalize(sequence3, 2);
        rnnPretrainingLabelMapper.mapLabels(sequence3, inputs, 2);
        rnnPretrainingLabelMapper.maskLabels(sequence3, mask, 2);

        assertEquals(inputs.toString(), expectedLabels);
        assertEquals(mask.toString(), expectedMask);
    }


}
