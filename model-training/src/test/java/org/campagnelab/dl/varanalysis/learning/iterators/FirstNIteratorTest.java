package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.junit.Test;
import org.nd4j.linalg.dataset.DataSet;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 7/14/16.
 */
public class FirstNIteratorTest {
    @Test
    public void testNIterator() throws IOException {
        int count=0;
        FirstNIterator it = new FirstNIterator(new BaseInformationIterator("sample_data/protobuf/genotypes_proto_test_mutated_randomized", 32,
                new SimpleFeatureCalculator(), new SimpleFeatureCalculator()), 10);
        while (it.hasNext()) {
            DataSet next = it.next();
            count+=next.numExamples();
        }
    assertEquals(32,count);
    }

    @Test
    public void testNIteratorBatch() throws IOException {
        int count=0;
        FirstNIterator it = new FirstNIterator(new BaseInformationIterator("sample_data/protobuf/genotypes_proto_test_mutated_randomized", 32,
                new SimpleFeatureCalculator(), new SimpleFeatureCalculator()), 10);
        while (it.hasNext()) {
            DataSet next = it.next(10);
            count+=next.numExamples();
        }
        assertEquals(10,count);
    }

}