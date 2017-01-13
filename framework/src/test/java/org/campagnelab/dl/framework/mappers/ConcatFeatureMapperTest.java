package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.function.Function;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 1/10/17.
 */
public class ConcatFeatureMapperTest {
    @Test
    public void concatFeatures() {
        String record = "012";
        String expectedLabels = "[1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00]";
        FeatureMapper<String>[] calculators = new FeatureMapper[record.length()];
        for (int i = 0; i < record.length(); i++) {
            calculators[i] = new OneHotBaseFeatureMapper<>(i, Function.identity(),
                    (r, idx) -> Character.getNumericValue(r.charAt(idx)), 3);
        }
        FeatureMapper<String> concatCalculator = new ConcatFeatureMapper<>(calculators);
        INDArray labels = Nd4j.zeros(1, concatCalculator.numberOfFeatures());
        concatCalculator.prepareToNormalize(record, 0);
        concatCalculator.mapFeatures(record, labels, 0);
        assertEquals(labels.toString(), expectedLabels);
    }
}
