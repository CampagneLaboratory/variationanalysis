package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Created by joc2080 on 11/28/16.
 */
public class RNNLabelMapperTest {
    @Test
    public void mapRNNLabelsOHBLMDelegates() {
        String sequence1 = "ATCGNJ";
        String sequence2 = "ATCG";
        String sequence3 = "AT";
        String expectedLabels =
                "[[[1.00, 0.00, 0.00, 0.00, 0.00, 1.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 1.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 1.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 1.00, 1.00, 0.00, 0.00]],\n" +
                "\n" +
                " [[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
                "  [0.00, 0.00, 0.00, 0.00, 0.00, 0.00]]]";
        String expectedMask =
                "[[1.00, 1.00, 1.00, 1.00, 1.00, 1.00],\n" +
                " [1.00, 1.00, 1.00, 1.00, 0.00, 0.00],\n" +
                " [1.00, 1.00, 0.00, 0.00, 0.00, 0.00]]";
        RNNLabelMapper<String> rnnLabelMapper = new RNNLabelMapper<>(6, 3,
                r -> Arrays.stream(r.split("")).mapToInt(this::testConvert).toArray(), String::length);

        INDArray inputs = Nd4j.zeros(3, 3, 6);
        INDArray mask = Nd4j.zeros(3, 6);

        rnnLabelMapper.prepareToNormalize(sequence1, 0);
        rnnLabelMapper.mapLabels(sequence1, inputs, 0);
        rnnLabelMapper.maskLabels(sequence1, mask, 0);

        INDArray sequence1Labels = inputs.getRow(0);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                int labelIndex = i * 3 + j;
                assertEquals(sequence1Labels.getFloat(j, i), rnnLabelMapper.produceLabel(sequence1, labelIndex), 1e-9);
            }
        }

        rnnLabelMapper.prepareToNormalize(sequence2, 1);
        rnnLabelMapper.mapLabels(sequence2, inputs, 1);
        rnnLabelMapper.maskLabels(sequence2, mask, 1);

        INDArray sequence2Mask = mask.getRow(1);
        for (int i = 0; i < 6; i++) {
            assertEquals(sequence2Mask.getInt(i) == 1, rnnLabelMapper.isMasked(sequence2, i * 3));
        }

        rnnLabelMapper.prepareToNormalize(sequence3, 2);
        rnnLabelMapper.mapLabels(sequence3, inputs, 2);
        rnnLabelMapper.maskLabels(sequence3, mask, 2);

        assertEquals(inputs.toString(), expectedLabels);
        assertEquals(mask.toString(), expectedMask);
    }

    private int testConvert(String oneChar) {
        switch (oneChar) {
            case "A":
            case "J":
                return 0;
            case "T":
            case "N":
                return 1;
            case "C":
            case "G":
                return 2;
            default:
                throw new RuntimeException("Shouldn't reach default");
        }
    }
}
