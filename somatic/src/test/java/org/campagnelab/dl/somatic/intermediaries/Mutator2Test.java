package org.campagnelab.dl.somatic.intermediaries;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import java.io.StringWriter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * Test the mutator on some specific examples.
 * Created by fac2003 on 5/27/16.
 */
public class Mutator2Test {

    @Test
    public void mutateTest() throws Exception {
        int index = 0;
        for (String record : records) {
            SimulationStrategy strategy = new SimulationStrategyImpl(0,1,0.1,1,0.9);
            strategy.setSeed(1);

            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);
            BaseInformationRecords.BaseInformation mutated = strategy.mutate(true, builder.build(), builder.getSamples(0), builder.getSamples(1), null);
            StringWriter writer = new StringWriter();
            TextFormat.print(mutated, writer);
            assertEquals("failed for record " + index, expectedMutatedRecords[index], writer.getBuffer().toString());
            index++;
        }
    }

    String hetExample = "reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: false\n" +
            "referenceBase: \"A\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 5\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +

            "} ";


    String[] records = {"reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: false\n" +
            "referenceBase: \"A\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 11\n" +
            "    genotypeCountReverseStrand: 14\n" +
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
            "    genotypeCountForwardStrand: 7\n" +
            "    genotypeCountReverseStrand: 13\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 400\n" +
            "    genotypeCountReverseStrand: 400\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 13\n" +
            "    genotypeCountReverseStrand: 16\n" +
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
            "    genotypeCountForwardStrand: 1\n" +
            "    genotypeCountReverseStrand: 2\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 400\n" +
            "    genotypeCountReverseStrand: 400\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "} ",

            "reference_index: 18\n" +
                    "position: 17214616\n" +
                    "mutated: false\n" +
                    "referenceBase: \"A\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 11\n" +
                    "    genotypeCountReverseStrand: 14\n" +
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
                    "    genotypeCountForwardStrand: 7\n" +
                    "    genotypeCountReverseStrand: 13\n" +
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
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 13\n" +
                    "    genotypeCountReverseStrand: 16\n" +
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
                    "    genotypeCountForwardStrand: 1\n" +
                    "    genotypeCountReverseStrand: 2\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "} ",

    };
    String[] expectedMutatedRecords = {"reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: true\n" +
            "mutatedBase: \"G\"\n" +
            "referenceBase: \"A\"\n" +
            "frequencyOfMutation: 0.7020648\n" +
            "indexOfMutatedBase: 3\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 11\n" +
            "    genotypeCountReverseStrand: 14\n" +
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
            "    genotypeCountForwardStrand: 7\n" +
            "    genotypeCountReverseStrand: 13\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 400\n" +
            "    genotypeCountReverseStrand: 400\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "  isTumor: false\n" +
            "}\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 13\n" +
            "    genotypeCountReverseStrand: 16\n" +
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
            "    genotypeCountForwardStrand: 130\n" +
            "    genotypeCountReverseStrand: 147\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"ABC\"\n" +
            "    toSequence: \"-\"\n" +
            "    genotypeCountForwardStrand: 271\n" +
            "    genotypeCountReverseStrand: 255\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "  isTumor: true\n" +
            "  formattedCounts: \"mutated (G) sample counts A=29 T=0 C=0 G=277 N=0 n/a indels:[526]\"\n" +
            "}\n"
            ,
            "reference_index: 18\n" +
                    "position: 17214616\n" +
                    "mutated: false\n" +   // not mutated because genotype indices differ between the two samples
                    "referenceBase: \"A\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 11\n" +
                    "    genotypeCountReverseStrand: 14\n" +
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
                    "    genotypeCountForwardStrand: 7\n" +
                    "    genotypeCountReverseStrand: 13\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  isTumor: false\n" +
                    "}\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 13\n" +
                    "    genotypeCountReverseStrand: 16\n" +
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
                    "    genotypeCountForwardStrand: 1\n" +
                    "    genotypeCountReverseStrand: 2\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  isTumor: true\n" +
                    "}\n"
    };

}