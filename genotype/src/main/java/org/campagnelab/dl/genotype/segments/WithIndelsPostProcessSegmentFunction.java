package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;


public class WithIndelsPostProcessSegmentFunction implements PostProcessSegmentFunction {
    @Override
    public Segment apply(Segment segment) {


        Iterable<BaseInformationRecords.BaseInformation> records = segment.getAllRecords();
        records.forEach((BaseInformationRecords.BaseInformation record) -> {
            int longestIndelLength = 0;
            String longestReference=record.getReferenceBase();
            for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                    if (count.getIsIndel()) {
                        final String fromSequence = count.getFromSequence();
                        longestIndelLength = Math.max(longestIndelLength, fromSequence.length());
                        longestIndelLength = Math.max(longestIndelLength, count.getToSequence().length());
                        if (fromSequence.length()>longestReference.length()) {
                            longestReference=fromSequence;
                        }
                    }
                }

            }
            if (longestIndelLength > 1) {
                for (int offset = 0; offset < longestIndelLength; offset++) {
                    BaseInformationRecords.BaseInformation.Builder copy = record.toBuilder();
                    //    System.out.printf("record position: %d %n",record.getPosition());
                    copy = segment.adjustCounts(copy, offset,longestReference);
                    segment.insertAfter(record, copy);
                }
                segment.hideRecord(record);
            }
        });


        return segment;
    }


}
