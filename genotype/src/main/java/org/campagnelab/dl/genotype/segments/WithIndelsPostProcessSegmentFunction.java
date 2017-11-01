package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;


public class WithIndelsPostProcessSegmentFunction implements java.util.function.Function<Segment, Segment> {
    @Override
    public Segment apply(Segment segment) {

        segment.populateTrueGenotypes();
        Iterable<BaseInformationRecords.BaseInformation> records = segment.getAllRecords();
        records.forEach(record -> {
            int longestIndelLength = 0;
            for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                    if (count.getIsIndel()) {
                        longestIndelLength = Math.max(longestIndelLength, count.getFromSequence().length());
                        longestIndelLength = Math.max(longestIndelLength, count.getToSequence().length());
                    }
                }

            }
            if (longestIndelLength > 1) {
                for (int offset = 0; offset < longestIndelLength; offset++) {
                    BaseInformationRecords.BaseInformation.Builder copy = record.toBuilder();
                    //    System.out.printf("record position: %d %n",record.getPosition());
                    copy = segment.recordList.adjustCounts(copy, offset);
                    segment.insertAfter(record, copy);
                }

                segment.recordList.hideRecord(record);
            }
        });


        return segment;
    }


}
