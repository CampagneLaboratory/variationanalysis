package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.genotype.segments.splitting.NoSplitStrategy;
import org.campagnelab.dl.genotype.segments.splitting.SingleCandidateIndelSegment;
import org.campagnelab.dl.genotype.segments.splitting.SingleCandidateIndelSplitStrategy;
import org.campagnelab.dl.genotype.segments.splitting.SplitStrategy;
import org.campagnelab.dl.genotype.tools.SBIToSSIConverterArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;

import static org.junit.Assert.assertEquals;

/**
 * Created by mas2182 on 10/30/17.
 */
@RunWith(JUnit4.class)
public class SplitStrategyTest {


    String expectedSnps ="ref=A\ttrueGenotype=A\tcounts= A=22  (from: A)\n"+
            "ref=A\ttrueGenotype=A/-\tcounts= A=32  -=33  (from: A)\n"+
            "ref=A\ttrueGenotype=A/T\tcounts= A=21  T=22  (from: A)\n"+
            "ref=A\ttrueGenotype=A/T\tcounts= A=21  T=22  (from: A)\n"+
            "ref=-\ttrueGenotype=-\tcounts= -=21  -=22  (from: -)\n"+
            "ref=C\ttrueGenotype=C\tcounts= C=15  (from: C)\n";

    SegmentHelper helper;

    @Before
    public void buildSegmentHelper() {
        Function<Segment, Segment> function = segment -> segment;
        SBIToSSIConverterArguments args = new SBIToSSIConverterArguments();
        args.mapFeatures = false;
        args.mapLabels = false;
        FillInFeaturesFunction fillInFeatures = new MyFillInFeaturesFunction(null, null, args);


        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segmentInfoo -> {
            //assertEquals(expectedSnps, Segment.showGenotypes(segmentInfoo));
        };
        helper = new SegmentHelper(function, fillInFeatures, segmentConsumer, new NoSplitStrategy(),
                true);

    }

    private BaseInformationRecords.BaseInformation makeIndel(int position) {
        return SegmentHelperTest.makeRecord(0,  position, "A--A", "A/A=20+12", "A/-=10+23");
    }

    private BaseInformationRecords.BaseInformation makeSnp(int position) {
        return SegmentHelperTest.makeRecord(0, position, "A/A", "A/A=12+10");

    }

    private void buildSBI(String sequence) {
        int position = 1;
        for (char c : sequence.toCharArray()){
            switch (c) {
                case 'S':
                    if (position==1)
                        helper.newSegment(makeSnp(position));
                    else
                        helper.add(makeSnp(position));

                    break;
                case 'I':
                    if (position==1)
                        helper.newSegment(makeIndel(position));
                    else
                        helper.add(makeIndel(position));

                    break;
            }
          position++;
        }
    }

    @Test
    public void testNoStrategy() {
        buildSBI("SIS");
        this.printCurrentIndels();
        SplitStrategy strategy = new NoSplitStrategy();
        List<Segment> subsegments = strategy.apply(helper.getCurrentSegment());
        assertEquals("Invalid number of subsegments returned by NoSplitStrategy", 1, subsegments.size());
    }

    @Test
    public void testSingleCandidateIndelSSISS() {
        buildSBI("SSISS");
        this.printCurrentIndels();
        SplitStrategy strategy = new SingleCandidateIndelSplitStrategy(1,0,true);
        List<SingleCandidateIndelSegment> subsegments = strategy.apply(helper.getCurrentSegment());
        assertEquals("Invalid number of subsegments returned by SingleCandidateIndelSplitStrategy", 1, subsegments.size());
        assertEquals("Invalid limits for the subsegment", "2-4", String.format("%d-%d",subsegments.get(0).getFirstPosition(),
                subsegments.get(0).getLastPosition()));
        assertEquals("Invalid indel position in the subsegment", 3, subsegments.get(0).getIndelPosition());
        helper.close();
    }

    @Test
    public void testSingleCandidateIndelSSIISS() {
        buildSBI("SSIISS");
        this.printCurrentIndels();
        SplitStrategy strategy = new SingleCandidateIndelSplitStrategy(1,0,true);
        List<SingleCandidateIndelSegment> subsegments = strategy.apply(helper.getCurrentSegment());
        assertEquals("Invalid number of subsegments returned by SingleCandidateIndelSplitStrategy", 1, subsegments.size());
        assertEquals("Invalid limits for the subsegment", "2-5", String.format("%d-%d",subsegments.get(0).getFirstPosition(),
                subsegments.get(0).getLastPosition()));
        assertEquals("Invalid indel position in the subsegment", 4, subsegments.get(0).getIndelPosition());
        helper.close();
    }

    @Test
    public void testSingleCandidateIndelSISISI() {
        buildSBI("SISISI");
        this.printCurrentIndels();
        SplitStrategy strategy = new SingleCandidateIndelSplitStrategy(1,0,true);
        List<SingleCandidateIndelSegment> subsegments = strategy.apply(helper.getCurrentSegment());
        assertEquals("Invalid number of subsegments returned by SingleCandidateIndelSplitStrategy", 3, subsegments.size());
        //first subsegment
        assertEquals("Invalid limits for the subsegment", "1-3", String.format("%d-%d",subsegments.get(0).getFirstPosition(),
                subsegments.get(0).getLastPosition()));
        assertEquals("Invalid indel position in the subsegment", 2, subsegments.get(0).getIndelPosition());
        //second subsegment
        assertEquals("Invalid limits for the subsegment", "3-5", String.format("%d-%d",subsegments.get(1).getFirstPosition(),
                subsegments.get(1).getLastPosition()));
        assertEquals("Invalid indel position in the subsegment", 4, subsegments.get(1).getIndelPosition());
        // third subsegment
        assertEquals("Invalid limits for the subsegment", "5-6", String.format("%d-%d",subsegments.get(2).getFirstPosition(),
                subsegments.get(2).getLastPosition()));
        assertEquals("Invalid indel position in the subsegment", 6, subsegments.get(2).getIndelPosition());
        helper.close();
    }

    private void printCurrentIndels() {
        Segment segment = helper.getCurrentSegment();
        Iterable<BaseInformationRecords.BaseInformation> it = segment.getAllRecords();
        it.forEach( record -> {
            System.out.println("Has candidate indel? " + SegmentUtil.hasCandidateIndel(record,0));
            System.out.println("Has true indel? " + SegmentUtil.hasTrueIndel(record));
        });
    }
}
