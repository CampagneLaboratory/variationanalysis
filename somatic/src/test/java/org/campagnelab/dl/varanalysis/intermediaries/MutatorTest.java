package org.campagnelab.dl.varanalysis.intermediaries;

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
public class MutatorTest {

    @Test
    public void mutateTest() throws Exception {
        int index = 0;
        for (String record : records) {
            Mutator m = new Mutator();
            m.setSeed(1);
            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);
            BaseInformationRecords.BaseInformation result = m.mutate(builder);
            StringWriter writer = new StringWriter();
            TextFormat.print(result, writer);
            assertEquals(expectedMutatedRecords[index], writer.getBuffer().toString());
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

    // disabled, the mutation is valid.
    public void testDontMutateAcrossHet() throws Exception {
        int index = 0;
        Mutator m = new Mutator();
        m.setSeed(1);
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(hetExample, builder);
        int MAX_TRIALS = 10000;
        for (int i = 0; i < MAX_TRIALS; i++) {
            BaseInformationRecords.BaseInformation.Builder germline = builder.clone();
            BaseInformationRecords.BaseInformation result = m.mutate(builder);
            BaseInformationRecords.SampleInfo germlineSample = germline.getSamples(0);
            int countExistingBase0= germlineSample.getCounts(0).getGenotypeCountForwardStrand()+ germlineSample.getCounts(0).getGenotypeCountReverseStrand();
            BaseInformationRecords.SampleInfo somaticSample = result.getSamples(1);
            int countExistingBase1= somaticSample.getCounts(1).getGenotypeCountForwardStrand()+ somaticSample.getCounts(1).getGenotypeCountReverseStrand();
            if (countExistingBase1>countExistingBase0) {
                System.out.println("germline: "+germline);
                System.out.println("somatic:  "+result);
                fail("Do not mutate bases towards pre-existing het base ");
            }else{
                // OK
            }
        }

    }

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
            "  isTumor: true\n"+
            "}\n"
            ,
            "reference_index: 18\n" +
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
                    "  isTumor: false\n" +
                    "}\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 7\n" +
                    "    genotypeCountReverseStrand: 10\n" +
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
                    "    genotypeCountReverseStrand: 8\n" +
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