package org.campagnelab.dl.somatic.learning.iterators;

import org.campagnelab.dl.somatic.mappers.IndelFeatures;
import org.campagnelab.dl.somatic.mappers.SimpleFeatureCalculator;
import org.junit.Test;
import org.nd4j.linalg.dataset.api.DataSet;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 11/2/16.
 */
public class CachedConcatIteratorTest {
    @Test
    public void test() throws IOException {
        // write a mapped features file:
        List<NamedDataSetIterator> iterators = new ArrayList<>();
        iterators.add(new BaseInformationIterator("sample_data/protobuf/genotypes_proto_mutated_randomized.sbi", 32, new IndelFeatures(), new SimpleFeatureCalculator()));
        CachedConcatIterator it = new CachedConcatIterator(iterators, 32, new IndelFeatures(), Integer.MAX_VALUE);
        long numRecords = 0;
        while (it.hasNext()) {
            DataSet ds = it.next();
            numRecords += ds.numExamples();
        }
        assertEquals(291884,numRecords);
        assertTrue("Cache file must exist",new File("sample_data/protobuf/genotypes_proto_mutated_randomized.cf").exists());
        assertTrue("Cache file properties must exist",new File("sample_data/protobuf/genotypes_proto_mutated_randomized.cfp").exists());

    }
}