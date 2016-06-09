package org.campagnelab.dl.varanalysis.learning.iterators;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.learning.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 5/27/16.
 */
public class SimpleFeatureCalculatorTest {
    @Test
    public void mapFeatures() throws Exception {
int index=0;

        for (String record : records) {
           SimpleFeatureCalculator calculator=new SimpleFeatureCalculator();

            INDArray inputs = Nd4j.zeros(1, calculator.numberOfFeatures());

            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);

            calculator.prepareToNormalize(builder.build(),0);
            calculator.mapFeatures(builder.build(),inputs,0);


            assertEquals(expectedFeatures[index], inputs.toString());
            index++;
        }
    }

    String[] records = {"reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: true\n" +
            "mutatedBase: \"C\"\n"+
            "indexOfMutatedBase: 2\n"+
            "referenceBase: \"A\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  isTumor: true\n"+
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "} "};
    String[] expectedFeatures = {"[0.12, 0.12, 0.00, 0.00, 0.00, 0.00, 0.12, 0.12, 0.00, 0.00, 0.12, 0.12, 0.00, 0.00, 0.06, 0.06, 0.06, 0.06, 0.00, 0.00, 0.12, 0.12, 0.12, 0.12, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12, 0.12, 0.06, 0.06, 0.06, 0.06, 0.00, 0.00, 0.00, 0.00]"};
}