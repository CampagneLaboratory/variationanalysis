package org.campagnelab.dl.varanalysis.tools;

import org.campagnelab.dl.varanalysis.learning.iterators.MappedFeaturesIterator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.dataset.api.DataSet;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 11/2/16.
 */
public class MapFeaturesTest {
    @Test
    public void test() {
        // write a mapped features file:
        MapFeatures tool = new MapFeatures();
        tool.arguments=new MapFeaturesArguments();
        tool.arguments.trainingSets.add("sample_data/protobuf/genotypes_proto_mutated_randomized.sbi");
        tool.arguments.miniBatchSize=32;
        tool.arguments.featureMapperClassname="org.campagnelab.dl.model.utils.mappers.IndelFeatures";
        tool.arguments.outputBasename="test-results/mapped-features";
        tool.execute();

        // now iterate through it:

        MappedFeaturesIterator it=new MappedFeaturesIterator(tool.arguments.outputBasename);
        long numRecords=0;
        while (it.hasNext()) {
            DataSet ds=it.next();
            numRecords+=ds.numExamples();
        }
        // Check we get back the same number of records written:
        assertEquals(tool.getNumRecordsWritten(),numRecords);
    }
}