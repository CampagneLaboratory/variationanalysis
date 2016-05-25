package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.varanalysis.learning.mappers.FeatureMapperV2;
import org.junit.Test;
import org.nd4j.linalg.dataset.DataSet;

import java.io.IOException;

import static org.junit.Assert.*;

/**
 * Run simple tests on the iterator.
 * Created by fac2003 on 5/21/16.
 */
public class BaseInformationIteratorTest {
    @Test
    public void testIterator() throws IOException{

        BaseInformationIterator trainIter = new BaseInformationIterator("test-data/genotypes_mutated_protofbuf.parquet",
                2, new SimpleFeatureCalculator(),new SimpleFeatureCalculator());
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        System.out.println(dataset );
    }

    @Test
    public void testIteratorV2() throws IOException{

        BaseInformationIterator trainIter = new BaseInformationIterator("test-data/genotypes_mutated_protofbuf.parquet",
                2, new FeatureMapperV2(), new SimpleFeatureCalculator());
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        System.out.println(dataset );
    }
}