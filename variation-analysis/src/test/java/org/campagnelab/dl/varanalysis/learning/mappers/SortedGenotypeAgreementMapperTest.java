package org.campagnelab.dl.varanalysis.learning.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 6/10/16.
 */
public class SortedGenotypeAgreementMapperTest {
    //update this test for new sort
    //Test
    public void mapFeatures() throws Exception {
        int index = 0;

        for (String record : records) {
            SortedGenotypeAgreementMapper calculator = new SortedGenotypeAgreementMapper();

            INDArray inputs = Nd4j.zeros(1, calculator.numberOfFeatures());

            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);

            calculator.prepareToNormalize(builder.build(), 0);
            calculator.mapFeatures(builder.build(), inputs, 0);


            assertEquals(String.format("test %d is expected to match", index), expectedFeatures[index], inputs.toString());
            index++;
        }
    }

    String[] records = {"reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: true\n" +
            "mutatedBase: \"C\"\n" +
            "indexOfMutatedBase: 2\n" +
            "referenceBase: \"A\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 10\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 9\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 8\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 7\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 6\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  isTumor: true\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 9\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 8\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 7\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 6\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "} " /* end of record 0 */
            ,
            "reference_index: 18\n" +
                    "position: 17214616\n" +
                    "mutated: true\n" +
                    "mutatedBase: \"C\"\n" +
                    "indexOfMutatedBase: 2\n" +
                    "referenceBase: \"A\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 150\n" +
                    "    genotypeCountReverseStrand: 150\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 120\n" +
                    "    genotypeCountReverseStrand: 120\n" +
                    "  }\n" +
                    "}\n" +
                    "samples {\n" +
                    "  isTumor: true\n" +
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
                    "    genotypeCountForwardStrand: 110\n" +
                    "    genotypeCountReverseStrand: 110\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +

                    "} " /* end of record 1 */
            ,
            "reference_index: 18\n" +
                    "position: 17214616\n" +
                    "mutated: true\n" +
                    "mutatedBase: \"C\"\n" +
                    "indexOfMutatedBase: 2\n" +
                    "referenceBase: \"A\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 10\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 9\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "}\n" +
                    "samples {\n" +
                    "  isTumor: true\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 9\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 10\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "} "
    };
    String[] expectedFeatures = {"[1.00, 1.00, 1.00, 1.00, 1.00]", "[0.00, 0.00, 0.00, 1.00, 1.00]", "[0.00, 0.00, 1.00, 1.00, 1.00]"};

}