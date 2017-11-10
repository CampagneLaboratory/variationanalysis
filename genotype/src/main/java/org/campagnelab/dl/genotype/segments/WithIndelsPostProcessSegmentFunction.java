package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.ArrayList;
import java.util.List;


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
            BaseInformationRecords.BaseInformation.Builder recordCopyWithIndices = record.toBuilder();
            if (longestIndelLength > 1) {
                List<Integer> indices = new IntArrayList();
                int i = 0;
                for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                    List<Integer> sampleIndices = new IntArrayList();
                    for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                        if (count.hasIsCalled() && count.getIsCalled() && count.hasGobyGenotypeIndex()) {
                            indices.add(count.getGobyGenotypeIndex());
                            sampleIndices.add(count.getGobyGenotypeIndex());
                        }
                    }
                    BaseInformationRecords.SampleInfo sampleWithIndices  = sample.toBuilder()
                            .addAllIndices(sampleIndices)
                            .build();
                    recordCopyWithIndices.setSamples(i, sampleWithIndices);
                    i++;
                }
                record = recordCopyWithIndices.build();
                for (int offset = 0; offset < longestIndelLength; offset++) {
                    BaseInformationRecords.BaseInformation.Builder copy = record.toBuilder();
                    //    System.out.printf("record position: %d %n",record.getPosition());
                    copy.addAllIndices(indices);
                    copy = segment.adjustCounts(copy, offset,longestReference);
                    segment.insertAfter(record, copy);
                }
                segment.hideRecord(record);
            }
        });


        return segment;
    }


}
