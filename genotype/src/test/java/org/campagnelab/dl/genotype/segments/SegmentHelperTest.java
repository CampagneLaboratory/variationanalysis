package org.campagnelab.dl.genotype.segments;


import org.campagnelab.dl.genotype.segments.splitting.NoSplitStrategy;
import org.campagnelab.dl.genotype.tools.SBIToSSIConverterArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.junit.Test;

import java.util.function.Consumer;
import java.util.function.Function;

import static org.junit.Assert.assertEquals;

public class SegmentHelperTest {
    String expected1 = "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
            "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
            "ref=C\ttrueGenotype=C\tcounts= C=15  (from: C)\n" ;


    @Test
    public void testHelper1() {
        Function<Segment, Segment> function = segment -> segment;
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expected1, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex=0;
        int position=0;
        helper.add(makeRecord(refIndex,position,"A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex,position+1,"A/T", "A/A=20+12", "A/T=10+23"));
    //    helper.add(makeRecord(refIndex,position+2,"AA/A--A", "AA/A--A=2+5","AA/AA=2+5"));
        helper.add(makeRecord(refIndex,position+3,"C/C", "C/C=10+5"));
        helper.close();

    }

    // format of count creation instruction is from/to=10+12
    private BaseInformationRecords.BaseInformation makeRecord(int refIndex, int position, String genotype, String... countCreations) {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setTrueGenotype(genotype);
        builder.setReferenceIndex(refIndex);
        builder.setPosition(position);

        BaseInformationRecords.SampleInfo.Builder sample = BaseInformationRecords.SampleInfo.newBuilder();
        for (String countCreationInstruction : countCreations) {
            BaseInformationRecords.CountInfo.Builder countBuilder = BaseInformationRecords.CountInfo.newBuilder();
            String tokens[] = countCreationInstruction.split("[/=+]");
            assert tokens.length == 4 :
                    "count creation instruction must have four arguments: ref/to=forward+reverse, was "+countCreationInstruction;
            countBuilder.setFromSequence(tokens[0]);
            builder.setReferenceBase(Character.toString(tokens[0].charAt(0)));
            countBuilder.setToSequence(tokens[1]);
            countBuilder.setMatchesReference(tokens[0].equals(tokens[1]));
            countBuilder.setGenotypeCountForwardStrand(Integer.parseInt(tokens[2]));
            countBuilder.setGenotypeCountReverseStrand(Integer.parseInt(tokens[3]));

            sample.addCounts(countBuilder);
        }
        sample.setFormattedCounts(FormatterCountHelper.format(sample));
        builder.addSamples(sample);
        return builder.build();
    }
}