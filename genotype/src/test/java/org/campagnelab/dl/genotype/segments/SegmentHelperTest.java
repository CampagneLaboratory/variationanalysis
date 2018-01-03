package org.campagnelab.dl.genotype.segments;


import org.campagnelab.dl.genotype.segments.splitting.NoSplitStrategy;
import org.campagnelab.dl.genotype.segments.splitting.SingleCandidateIndelSplitStrategy;
import org.campagnelab.dl.genotype.tools.SBIToSSIConverterArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.junit.Test;

import java.util.function.Consumer;
import java.util.function.Function;

import static org.junit.Assert.assertEquals;

public class SegmentHelperTest {
    String expectedSnps =
            "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
            "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
            "ref=C\ttrueGenotype=C\tcounts= C=15  (from: C)\n";

    String expectedHetInsertionCC =
            "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
                    "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
                    "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n" +
                    "ref=-\ttrueGenotype=-/C\tcounts= -=7  C=7  (from: -)\n" +
                    "ref=-\ttrueGenotype=-/C\tcounts= -=7  C=7  (from: -)\n" +
                    "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n";

    String expectedHomInsertionCC =
            "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
        "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
        "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n" +
        "ref=-\ttrueGenotype=C\tcounts= C=14  (from: -)\n" +
        "ref=-\ttrueGenotype=C\tcounts= C=14  (from: -)\n" +
        "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n";

    String expectedHetDeletionTT =
            "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n" +
            "ref=T\ttrueGenotype=-/T\tcounts= -=7  T=7  (from: T)\n" +
            "ref=T\ttrueGenotype=-/T\tcounts= -=7  T=7  (from: T)\n" +
            "ref=A\ttrueGenotype=A\tcounts= A=14  (from: A)\n";

    String expectedHomDeletionTT =
            "ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n" +
           "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A)\n" +
           "ref=A\ttrueGenotype=A\tcounts= A=7  (from: A)\n" +
           "ref=T\ttrueGenotype=-\tcounts= -=7  (from: T)\n" +
           "ref=T\ttrueGenotype=-\tcounts= -=7  (from: T)\n" +
           "ref=A\ttrueGenotype=A\tcounts= A=7  (from: A)\n";

    String expectedInsertionWithSNP =
            "ref=A\ttrueGenotype=A\tcounts= A=65  (from: A) position=0 offset=0 \n" +
            "ref=-\ttrueGenotype=-/T\tcounts= -=32  T=33  (from: -) position=0 offset=1 \n" +
            "ref=-\ttrueGenotype=-/T\tcounts= -=32  T=33  (from: -) position=0 offset=2 \n" +
            "ref=A\ttrueGenotype=A/T\tcounts= A=32  T=33  (from: A) position=1 offset=3 \n";
    String expectedInsertionWithSNP2 = expectedInsertionWithSNP;

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
    // het insertion
    public void testHelperInsertion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedHetInsertionCC, Segment.showGenotypes(segmentInfoo));
        };
        SingleCandidateIndelSplitStrategy strategy = new SingleCandidateIndelSplitStrategy(2,0,false);
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, strategy,
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));
        helper.add(makeRecord(refIndex, position + 2, "A--A/ACCA", "A--A/ACCA=2+5", "A--A/A--A=2+5"));
        //  helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
        helper.close();

    }

    @Test
    // hom insertion
    public void testHelperHomInsertion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedHomInsertionCC, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));
        helper.add(makeRecord(refIndex, position + 2, "ACCA/ACCA", "A--A/ACCA=2+5", "A--A/ACCA=2+5"));
        //  helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
        helper.close();

    }

    @Test
    // homozygous deletion
    public void testHelperHomDeletion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedHomDeletionTT, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));
        helper.add(makeRecord(refIndex, position + 2, "A--A/A--A", "ATTA/A--A=2+5"));
        //  helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
        helper.close();

    }

    @Test
    // het deletion
    public void testHelperHetDeletion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);

        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedHetDeletionTT, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
     //   helper.add(makeRecord(refIndex, position, "A/A", "A/A=12+10"));
      //  helper.add(makeRecord(refIndex, position + 1, "A/T", "A/A=20+12", "A/T=10+23"));
        helper.add(makeRecord(refIndex, position + 2, "A--A/ATTA", "ATTA/A--A=2+5","ATTA/ATTA=2+5"));
        //  helper.add(makeRecord(refIndex, position + 3, "C/C", "C/C=10+5"));
        helper.close();

    }

    // insertion w/ SNP
    public void testInsertionWithSNP() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);

        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedInsertionWithSNP, Segment.showGenotypes(segmentInfoo,true));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A--A/ATTT","A--A/A--A=20+12", "A--A/ATTT=10+23"));
        helper.close();

    }

    // insertion w/ SNP 2
    public void testInsertionWithSNP2() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);

        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedInsertionWithSNP2, Segment.showGenotypes(segmentInfoo,true));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                false);
        int refIndex = 0;
        int position = 0;
        helper.add(makeRecord(refIndex, position, "A--A/ATTA","A--A/A--A=20+12", "A--A/ATTA=10+23"));
        helper.add(makeRecord(refIndex, position + 1, "A/T","A/A=20+12", "A/T=10+23"));
        helper.close();

    }

    public void testBigInsertion() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);
        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            String segmentGenotypes = Segment.showGenotypes(segmentInfoo,true);
            assertEquals(expectedInsertionWithSNP2, segmentGenotypes);
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new SingleCandidateIndelSplitStrategy(10, 0, false),
                false);
        int refIndex = 0;
        int position = 112428313;
        helper.add(makeRecord(refIndex, position, "GTTTTTTTTTT/G-TTTTTTTTT","GTTTTTTTTTT/GTTTTTTTTTT=100+163", "GTTTTTTTTTT/G-TTTTTTTTT=200+137"));
        helper.close();
    }

    public void testBigInsertion2() {
        Function<Segment, Segment> function = new WithIndelsPostProcessSegmentFunction();
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);
        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            String segmentGenotypes = Segment.showGenotypes(segmentInfoo,true);
            assertEquals(expectedInsertionWithSNP2, segmentGenotypes);
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new SingleCandidateIndelSplitStrategy(10, 0, false),
                false);
        int refIndex = 0;
        int position = 112428313;
        helper.add(makeRecord(refIndex, position, "GTTTTTTTTTT/G-TTTTTTTTT",
                "GTTTTTTTTTT/G-TTTTTTTTT=183+111",
                "G/G=0+0",
                "GTTTTTTTTTT/GTTTTTTTTTT=106+144",
                "G/C=0+0",
                "G/N=0+0",
                "G/T=0+0",
                "G/A=0+0"));
        helper.close();
    }

    // format of count creation instruction is from/to=10+12
    protected static BaseInformationRecords.BaseInformation makeRecord(int refIndex, int position, String genotype, String... countCreations) {
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