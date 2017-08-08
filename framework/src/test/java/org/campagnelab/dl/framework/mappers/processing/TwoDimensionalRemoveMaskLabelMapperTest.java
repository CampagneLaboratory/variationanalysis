package org.campagnelab.dl.framework.mappers.processing;

import org.campagnelab.dl.framework.mappers.TwoDimensionalConcatLabelMapper;
import org.campagnelab.dl.framework.mappers.RNNLabelMapper;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 12/13/16.
 */
public class TwoDimensionalRemoveMaskLabelMapperTest {
    @Test
    public void testLabelFilterConcatRNNLabelMapperOHBMDelegates() {
        String sequence1 = "ATC";
        String sequence2 = "GN";
        String sequence3 = "J";
        String expectedLabels =
                "[[[1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00],  \n" +
                "  [0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00],  \n" +
                "  [0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00]],  \n" +
                "\n" +
                " [[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  1.00,  0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00],  \n" +
                "  [1.00,  0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00]],  \n" +
                "\n" +
                " [[1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask =
                "[[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00],  \n" +
                " [1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00],  \n" +
                " [1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]";
        RNNLabelMapper<String> rnnLabelMapper1 = new RNNLabelMapper<>(3, 3,
                r -> Arrays.stream(r.split("")).mapToInt(this::testConvert).toArray(), String::length);
        RNNLabelMapper<String> rnnLabelMapper2 = new RNNLabelMapper<>(3, 3,
                r -> Arrays.stream(r.split("")).mapToInt(this::testConvert).toArray(), String::length);
        RNNLabelMapper<String> rnnLabelMapper3 = new RNNLabelMapper<>(3, 3,
                r -> Arrays.stream(r.split("")).mapToInt(this::testConvert).toArray(), String::length);
        TwoDimensionalConcatLabelMapper<String> twoDConcatLMapper = new TwoDimensionalConcatLabelMapper<>(rnnLabelMapper1, rnnLabelMapper2, rnnLabelMapper3);
        TwoDimensionalRemoveMaskLabelMapper<String> twoDFilterLMapper = new TwoDimensionalRemoveMaskLabelMapper<>(twoDConcatLMapper);

        INDArray inputs = Nd4j.zeros(3, 3, 9);
        INDArray mask = Nd4j.zeros(3, 9);

        twoDFilterLMapper.prepareToNormalize(sequence1, 0);
        twoDFilterLMapper.mapLabels(sequence1, inputs, 0);
        twoDFilterLMapper.maskLabels(sequence1, mask, 0);

        INDArray sequence1Labels = inputs.getRow(0);
        for (int i = 0; i < 9; i++) {
            for (int j = 0; j < 3; j++) {
                int labelIndex = i * 3 + j;
                assertEquals(sequence1Labels.getFloat(j, i), twoDFilterLMapper.produceLabel(sequence1, labelIndex), 1e-9);
            }
        }

        twoDFilterLMapper.prepareToNormalize(sequence2, 1);
        twoDFilterLMapper.mapLabels(sequence2, inputs, 1);
        twoDFilterLMapper.maskLabels(sequence2, mask, 1);

        INDArray sequence2Mask = mask.getRow(1);
        for (int i = 0; i < 9; i++) {
            assertEquals(sequence2Mask.getInt(i) == 1, twoDFilterLMapper.isMasked(sequence2, i * 3));
        }

        twoDFilterLMapper.prepareToNormalize(sequence3, 2);
        twoDFilterLMapper.mapLabels(sequence3, inputs, 2);
        twoDFilterLMapper.maskLabels(sequence3, mask, 2);

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
