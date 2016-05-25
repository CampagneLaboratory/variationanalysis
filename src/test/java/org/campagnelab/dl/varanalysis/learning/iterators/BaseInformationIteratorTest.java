package org.campagnelab.dl.varanalysis.learning.iterators;

import org.junit.Test;
import org.nd4j.linalg.dataset.DataSet;

import static org.junit.Assert.*;

/**
 * Run simple tests on the iterator.
 * Created by fac2003 on 5/21/16.
 */
public class BaseInformationIteratorTest {
    @Test
    public void testIterator() {

        BaseInformationIterator trainIter = new BaseInformationIterator("sample_data/genotypes_mutated.parquet",
                2, new SimpleFeatureCalculator(),new SimpleFeatureCalculator());
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        System.out.println(dataset );
    }
}