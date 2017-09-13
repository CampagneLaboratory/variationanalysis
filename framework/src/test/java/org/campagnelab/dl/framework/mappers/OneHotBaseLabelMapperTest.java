package org.campagnelab.dl.framework.mappers;

import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 11/28/2016
 */
public class OneHotBaseLabelMapperTest {
    @Test
    public void mapLabels() throws Exception {
        for (int i = 0; i < record.length(); i++) {
            LabelMapper<String> calculator = new OneHotBaseLabelMapper<>(i, record.length(),
                    r -> Arrays.stream(r.split("")).mapToInt(Integer::parseInt).toArray());
            INDArray labels = Nd4j.zeros(1, calculator.numberOfLabels());
            calculator.prepareToNormalize(record,  0);
            calculator.mapLabels(record, labels, 0);
            assertEquals(expectedLabels[i], labels.toString());
        }
    }

    private String record = "012345";
    private String[] expectedLabels = {
            "[1.00,  0.00,  0.00,  0.00,  0.00,  0.00]",
            "[0.00,  1.00,  0.00,  0.00,  0.00,  0.00]",
            "[0.00,  0.00,  1.00,  0.00,  0.00,  0.00]",
            "[0.00,  0.00,  0.00,  1.00,  0.00,  0.00]",
            "[0.00,  0.00,  0.00,  0.00,  1.00,  0.00]",
            "[0.00,  0.00,  0.00,  0.00,  0.00,  1.00]",
    };
}