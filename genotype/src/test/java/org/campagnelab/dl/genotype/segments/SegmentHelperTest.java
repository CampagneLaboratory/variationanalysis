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
    String expectedSnps = "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
            "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
            "ref=C\ttrueGenotype=C\tcounts= C=15  (from: C)\n";

    String expectedInsertionCC = "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
            "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
            "ref=A\ttrueGenotype=A\tcounts= A=7  A=7  (from: A)\n" +
            "ref=-\ttrueGenotype=-/C\tcounts= C=7  -=7  (from: -)\n" +
            "ref=-\ttrueGenotype=-/C\tcounts= C=7  -=7  (from: -)\n" +
            "ref=A\ttrueGenotype=A\tcounts= A=7  A=7  (from: A)\n";


    @Test
    public void testHelperSnps() {
        Function<Segment, Segment> function = segment -> segment;
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedSnps, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));

        helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
        helper.close();

    }

    @Test
    public void testHelperInsertion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedInsertionCC, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));
        helper.add(makeRecord(refIndex, position + 2, "A--A/ACCA", "A--A/ACCA=2+5", "A--A/A--A=2+5"));
        //  helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
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
                    "count creation instruction must have four arguments: ref/to=forward+reverse, was " + countCreationInstruction;
            final String from = tokens[0];
            countBuilder.setFromSequence(from);
            builder.setReferenceBase(Character.toString(from.charAt(0)));
            final String token = tokens[1];
            countBuilder.setToSequence(token);
            countBuilder.setMatchesReference(from.equals(token));
            countBuilder.setGenotypeCountForwardStrand(Integer.parseInt(tokens[2]));
            countBuilder.setGenotypeCountReverseStrand(Integer.parseInt(tokens[3]));
            if (from.contains("-") || token.contains("-")) {
                countBuilder.setIsIndel(true);
            }

            sample.addCounts(countBuilder);
        }
        sample.setFormattedCounts(FormatterCountHelper.format(sample));
        builder.addSamples(sample);
        return builder.build();
    }
}