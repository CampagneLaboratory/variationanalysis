package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.genotype.segments.splitting.NoSplitStrategy;
import org.campagnelab.dl.genotype.segments.splitting.SplitStrategy;
import org.campagnelab.dl.genotype.tools.SBIToSSIConverterArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.junit.Test;

import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;

import static org.junit.Assert.assertEquals;

/**
 * Created by mas2182 on 10/30/17.
 */
public class SplitStrategyTest {


    String expectedSnps ="ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n"+
            "ref=A\ttrueGenotype=A/-\tcounts= A=32  -=33  (from: A)\n"+
            "ref=A\ttrueGenotype=A/T\tcounts= A=21  T=22  (from: A)\n"+
            "ref=A\ttrueGenotype=A/T\tcounts= A=21  T=22  (from: A)\n"+
            "ref=-\ttrueGenotype=-\tcounts= -=21  -=22  (from: -)\n"+
            "ref=C\ttrueGenotype=C\tcounts= C=15  (from: C)\n";


    @Test
    public void testNoStrategy() {

        Function<Segment, Segment> function = segment -> segment;
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            assertEquals(expectedSnps, Segment.showGenotypes(segmentInfoo));
        };
        SegmentHelper helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                true);
        int refIndex = 0;
        int position = 0;
        helper.add(SegmentHelperTest.makeRecord(refIndex, position, "A/A", "A/A=12+10"));
        helper.add(SegmentHelperTest.makeRecord(refIndex, position + 1, "A/-", "A/A=20+12", "A/-=10+23"));
        helper.add(SegmentHelperTest.makeRecord(refIndex, position + 4, "A/T", "A/A=10+11", "A/T=11+11"));
        helper.add(SegmentHelperTest.makeRecord(refIndex, position + 14, "A/T", "A/A=10+11", "A/T=11+11"));
        helper.add(SegmentHelperTest.makeRecord(refIndex, position + 24, "-/-", "-/-=10+11", "-/-=11+11"));
        helper.add(SegmentHelperTest.makeRecord(refIndex, position + 30, "C/C", "C/C=10+5"));
        Segment segment = helper.getCurrentSegment();
        Iterable<BaseInformationRecords.BaseInformation> it = segment.getAllRecords();
        it.forEach( record -> {
            System.out.println("Has candidate indel? " + SegmentUtil.hasCandidateIndel(record,0));
            System.out.println("Has true indel? " + SegmentUtil.hasTrueIndel(record));
        });
        helper.close();
        helper.printStats();
        SplitStrategy strategy = new NoSplitStrategy();
        List<Segment> subsegments = strategy.apply(segment);
        assertEquals("Invalid number of subsegments returned by NoSplitStrategu", 2, subsegments.size());

        System.out.println("Helper done.");
    }
}
