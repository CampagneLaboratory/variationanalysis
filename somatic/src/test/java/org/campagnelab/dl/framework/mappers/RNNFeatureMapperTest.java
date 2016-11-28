package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

import static org.junit.Assert.assertEquals;


/**
 * Created by joc2080 on 11/28/16.
 */
public class RNNFeatureMapperTest {
    @Test
    public void mapRNNFeatures() {
        RNNFeatureMapper<String> rnnFeatureMapper = new RNNFeatureMapper<>(6, 6, Function.identity(), String::length);
        INDArray inputs = Nd4j.zeros(1, 6, 6);
        rnnFeatureMapper.mapFeatures(sequence, inputs, 0);
        assertEquals(inputs.toString(), expectedFeature);
    }
    String sequence = "ATCGNJ";
    String expectedFeature = "[[[1.00, 0.00, 0.00, 0.00, 0.00, 0.00],\n" +
            "  [0.00, 1.00, 0.00, 0.00, 0.00, 0.00],\n" +
            "  [0.00, 0.00, 1.00, 0.00, 0.00, 0.00],\n" +
            "  [0.00, 0.00, 0.00, 1.00, 0.00, 0.00],\n" +
            "  [0.00, 0.00, 0.00, 0.00, 1.00, 0.00],\n" +
            "  [0.00, 0.00, 0.00, 0.00, 0.00, 1.00]]]";
}
