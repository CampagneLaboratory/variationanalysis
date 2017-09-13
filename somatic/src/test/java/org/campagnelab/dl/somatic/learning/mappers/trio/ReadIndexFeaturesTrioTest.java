package org.campagnelab.dl.somatic.learning.mappers.trio;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.somatic.mappers.trio.ReadIndexFeaturesTrio;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.assertEquals;

/**
 * Created by rct66 on 7/29/16.
 */
public class ReadIndexFeaturesTrioTest {
    @Test
    public void mapFeatures() throws Exception {
        int index=0;

        for (String record : records) {
            ReadIndexFeaturesTrio calculator = new ReadIndexFeaturesTrio();

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
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 5\n" +
            "    isIndel: true\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 5\n" +
            "    isIndel: true\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 98\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 1\n" +
            "      frequency: 2\n" +
            "    }\n" +
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
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 999\n" +
            "    genotypeCountReverseStrand: 999\n" +
            "    isIndel: true\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 99\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 2\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "  }\n" +
            "} "};
    String[] expectedFeatures = {"[0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.33,  0.33,  0.33,  0.33,  0.33]"};
}