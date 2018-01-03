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
            if (longestIndelLength > 1) {
                List<Integer> indices = new IntArrayList();
                List<List<Integer>> sampleIndicesList = new ArrayList<>();
                for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                    List<Integer> sampleIndices = new IntArrayList();
                    int gobyGenotypeIndex = 0;
                    for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                        if (count.hasIsCalled() && count.getIsCalled()) {
                            indices.add(gobyGenotypeIndex);
                            sampleIndices.add(gobyGenotypeIndex);
                        }
                        gobyGenotypeIndex++;
                    }
                    sampleIndicesList.add(sampleIndices);
                }
                for (int offset = 0; offset < longestIndelLength; offset++) {
                    //    System.out.printf("record position: %d %n",record.getPosition());
                    BaseInformationRecords.BaseInformation.Builder copy = record.toBuilder();
                    copy.addAllCalledCountIndices(indices);
                    int sampleIndex = 0;
                    for (BaseInformationRecords.SampleInfo.Builder sampleInfo : copy.getSamplesBuilderList()) {
                        sampleInfo.addAllCalledCountIndices(sampleIndicesList.get(sampleIndex));
                        copy.setSamples(sampleIndex, sampleInfo);
                        sampleIndex++;
                    }
                    copy = segment.adjustCounts(copy, offset, longestReference);
                    segment.insertAfter(record, copy);
                }
                segment.hideRecord(record);
                segment.setIndicesAdded(true);
            }
        });


        return segment;
    }

}
