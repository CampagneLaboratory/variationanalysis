package org.campagnelab.dl.somatic.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import static junit.framework.TestCase.assertEquals;


/**
 * Created by fac2003 on 12/5/16.
 */
public class IsBaseMutatedMapperTest {
    @Test
    public void testLabelMapper() throws TextFormat.ParseException {
        IsBaseMutatedMapper mapper = new IsBaseMutatedMapper(2);
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
      TextFormat.getParser().merge(proto, builder);
        BaseInformationRecords.BaseInformation record =builder.build();
        mapper.prepareToNormalize(record, 0);
        assertEquals(0f, mapper.produceLabel(record, 0));
        assertEquals(0f, mapper.produceLabel(record, 1));
        assertEquals(0f, mapper.produceLabel(record, 2));
        assertEquals(1f, mapper.produceLabel(record, 3)); // T is sorted in third base position.
        assertEquals(0f, mapper.produceLabel(record, 4));

        builder.setMutated(false);
        builder.clearMutatedBase();
        record=builder.build();
        mapper.prepareToNormalize(record, 0);
        assertEquals(1f, mapper.produceLabel(record, 0));// the new record is not mutated.
        assertEquals(0f, mapper.produceLabel(record, 1));
        assertEquals(0f, mapper.produceLabel(record, 2));
        assertEquals(0f, mapper.produceLabel(record, 3));
        assertEquals(0f, mapper.produceLabel(record, 4));

    }

    String proto = "reference_index: 18\n" +
            "position: 18282401\n" +
            "mutated: true\n" +
            "mutatedBase: \"T\"\n" +
            "referenceBase: \"G\"\n" +
            "frequencyOfMutation: 0.23167634\n" +
            "indexOfMutatedBase: 1\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 1\n" +
            "    isIndel: false\n" +

            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 7\n" +
            "    genotypeCountReverseStrand: 6\n" +
            "    isIndel: false\n" +

            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"GAAAAAAAAAAAAAAAAAA\"\n" +
            "    genotypeCountForwardStrand: 2\n" +
            "    genotypeCountReverseStrand: 2\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "  isTumor: false\n" +
            "  formattedCounts: \"sample: 1 counts A=1 T=0 C=0 G=13 N=0 FB=0 indels={ [indel count=2 G -AAAAAAAAAAAAAAAAA/AAAAAAAAAAAAAAAAAA  18282401-18282419 filtered=false] }\\n\"\n" +
            "}\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 2\n" +
            "    isIndel: false\n" +
            "    qualityScoresForwardStrand {\n" +
            "      number: 67\n" +
            "      frequency: 1\n" +
            "    }\n" +

            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 1\n" +
            "    genotypeCountReverseStrand: 2\n" +
            "    isIndel: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 20\n" +
            "    genotypeCountReverseStrand: 14\n" +
            "    isIndel: false\n" +

            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"G-AAAAAAAAAAAAAAAAA\"\n" +
            "    toSequence: \"GAAAAAAAAAAAAAAAAAA\"\n" +
            "    genotypeCountForwardStrand: 8\n" +
            "    genotypeCountReverseStrand: 7\n" +
            "    isIndel: true\n" +
            "  }\n" +
            "  isTumor: true\n" +
            "  formattedCounts: \"mutated (T) sample counts A=2 T=3 C=0 G=34 N=0 FB=0 indels:[15]\"\n" +
            "}\n" +
            "reference_id: \"19\"\n" +
            "genomicSequenceContext: \"TCTGTCTCAAGAAAAAAAAAA\"\n" ;


}