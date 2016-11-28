package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;
import java.util.function.Function;

import static org.junit.Assert.assertEquals;

/**
 * Created by joc2080 on 11/28/16.
 */
public class RNNLabelMapperTest {
    @Test
    public void mapRNNLabels() {
        RNNLabelMapper<String> rnnFeatureMapper = new RNNLabelMapper<>(3, 6,
                r -> Arrays.stream(r.split("")).mapToInt(this::testConvert).toArray(), String::length);
        INDArray labels = Nd4j.zeros(1, 3, 6);
        rnnFeatureMapper.mapLabels(sequence, labels, 0);
        assertEquals(labels.toString(), expectedLabel);
    }
    String sequence = "ATCGNJ";
    String expectedLabel = "[[[1.00, 0.00, 0.00, 0.00, 0.00, 1.00],\n" +
            "  [0.00, 1.00, 0.00, 0.00, 1.00, 0.00],\n" +
            "  [0.00, 0.00, 1.00, 1.00, 0.00, 0.00]]]";
    public int testConvert(String oneChar) {
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
