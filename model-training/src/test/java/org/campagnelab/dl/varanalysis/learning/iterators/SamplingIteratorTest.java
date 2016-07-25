package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.junit.Test;
import org.nd4j.linalg.dataset.DataSet;

import java.io.IOException;

import static org.junit.Assert.*;


/**
 * Created by fac2003 on 7/21/16.
 */
public class SamplingIteratorTest {
    @Test
    public void sample() throws IOException {
        long seed = 012;
        int minibatchSize = 3;
        SamplingIterator trainIter = new SamplingIterator(
                new BaseInformationIterator("sample_data/protobuf/genotypes_proto_mutated_randomized.parquet",
                        minibatchSize, new SimpleFeatureCalculator(), new SimpleFeatureCalculator()),
                seed);
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(minibatchSize, dataset.numExamples());
        for (int i = 0; i < dataset.numExamples(); i++) {

            trainIter.setSamplingProbability(true, i, 0.5f);
        }
        trainIter.reset();

        assertTrue(trainIter.hasNext());
        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(minibatchSize, dataset.numExamples());
        System.out.println(dataset);

    }
}