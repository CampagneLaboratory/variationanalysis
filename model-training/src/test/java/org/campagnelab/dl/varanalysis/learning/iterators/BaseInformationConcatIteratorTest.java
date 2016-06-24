package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV2;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.junit.Test;
import org.nd4j.linalg.dataset.DataSet;

import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Run simple tests on the concat iterator.
 * Created by rt66 on 6/14/16.
 */
public class BaseInformationConcatIteratorTest {
    @Test
    public void testIterator() throws IOException{

        List<BaseInformationIterator> iterList = new ObjectArrayList<BaseInformationIterator>(2);
        iterList.add( new BaseInformationIterator("sample_data/protobuf/concat_genotypes_test_proto_mutated.parquet",
                100, new SimpleFeatureCalculator(),new SimpleFeatureCalculator()));
        iterList.add( new BaseInformationIterator("sample_data/protobuf/concat_genotypes_test_proto_mutated.parquet",
                100, new SimpleFeatureCalculator(),new SimpleFeatureCalculator()));
        BaseInformationConcatIterator trainIter = new BaseInformationConcatIterator(iterList, 3, new FeatureMapperV2(), new SimpleFeatureCalculator());
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(3,dataset.numExamples());
        System.out.println(dataset );

        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(3,dataset.numExamples());
        System.out.println(dataset );

        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(2,dataset.numExamples());
        System.out.println(dataset );

        assert(!trainIter.hasNext());
    }

    @Test
    public void testIterator2() throws IOException{

        List<BaseInformationIterator> iterList = new ObjectArrayList<BaseInformationIterator>(2);
        iterList.add( new BaseInformationIterator("sample_data/protobuf/concat_genotypes_test_proto_mutated.parquet",
                100, new SimpleFeatureCalculator(),new SimpleFeatureCalculator()));
        iterList.add( new BaseInformationIterator("sample_data/protobuf/concat_genotypes_test_proto_mutated.parquet",
                100, new SimpleFeatureCalculator(),new SimpleFeatureCalculator()));
        BaseInformationConcatIterator trainIter = new BaseInformationConcatIterator(iterList, 2, new FeatureMapperV2(), new SimpleFeatureCalculator());
        assertTrue(trainIter.hasNext());
        DataSet dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(2,dataset.numExamples());
        System.out.println(dataset );

        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(2,dataset.numExamples());
        System.out.println(dataset );

        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(2,dataset.numExamples());
        System.out.println(dataset );

        dataset = trainIter.next();
        assertNotNull(dataset);
        assertEquals(2,dataset.numExamples());
        System.out.println(dataset );

        assert(!trainIter.hasNext());

    }

}