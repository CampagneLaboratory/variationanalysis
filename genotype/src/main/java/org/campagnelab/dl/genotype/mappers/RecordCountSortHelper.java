package org.campagnelab.dl.genotype.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Sort genotypes of a record by decreasing count.
 * Created by fac2003 on 12/15/16.
 */
public class RecordCountSortHelper {


    public BaseInformationRecords.BaseInformationOrBuilder sort(BaseInformationRecords.BaseInformationOrBuilder record) {
        int originalGenotypeIndex = 0;
        final List<BaseInformationRecords.CountInfo> countsList = record.getSamples(0).getCountsList();
        List<BaseInformationRecords.CountInfo> counts = new ArrayList<>();
        for (BaseInformationRecords.CountInfo count : countsList) {
            counts.add(BaseInformationRecords.CountInfo.newBuilder()
                    .mergeFrom(count)
                    .build());
        }

        Collections.sort(counts, (a, b) ->
                (b.getGenotypeCountForwardStrand() + b.getGenotypeCountReverseStrand()) -
                        (a.getGenotypeCountForwardStrand() + a.getGenotypeCountReverseStrand()));

        BaseInformationRecords.BaseInformation.Builder copyOfRecord = BaseInformationRecords.BaseInformation.newBuilder();
        copyOfRecord.setGenomicSequenceContext(record.getGenomicSequenceContext());
        copyOfRecord.setPosition(record.getPosition());
        copyOfRecord.setReferenceId(record.getReferenceId());
        copyOfRecord.setReferenceBase(record.getReferenceBase());
        copyOfRecord.setTrueGenotype(record.getTrueGenotype());
        copyOfRecord.setReferenceIndex(record.getReferenceIndex());
        final BaseInformationRecords.SampleInfo.Builder builder = record.getSamples(0).toBuilder();
        builder.clearCounts();
        builder.addAllCounts(counts);
        copyOfRecord.addSamples(builder);
        return copyOfRecord.build();
    }
    public BaseInformationRecords.BaseInformation sort(BaseInformationRecords.BaseInformation record) {
        final List<BaseInformationRecords.CountInfo> countsList = record.getSamples(0).getCountsList();
        List<BaseInformationRecords.CountInfo> counts = new ArrayList<>();
        for (BaseInformationRecords.CountInfo count : countsList) {
            int forwardCount = count.getIsIndel()
                    ? count.getGenotypeCountForwardStrand() - (count.getGenotypeCountForwardStrand() / 2)
                    : count.getGenotypeCountForwardStrand();
            int reverseCount = count.getIsIndel()
                    ? count.getGenotypeCountReverseStrand() / 2
                    : count.getGenotypeCountReverseStrand();
            counts.add(BaseInformationRecords.CountInfo.newBuilder()
                    .mergeFrom(count)
                    .setGenotypeCountForwardStrand(forwardCount)
                    .setGenotypeCountReverseStrand(reverseCount)
                    .build());
        }

        Collections.sort(counts, (a, b) ->
                (b.getGenotypeCountForwardStrand() + b.getGenotypeCountReverseStrand()) -
                        (a.getGenotypeCountForwardStrand() + a.getGenotypeCountReverseStrand()));


        final BaseInformationRecords.SampleInfo.Builder builder = record.getSamples(0).toBuilder();
        builder.clearCounts();
        builder.addAllCounts(counts);
        BaseInformationRecords.BaseInformation.Builder copyOfRecord = record.toBuilder();
        copyOfRecord.setSamples(0,builder);
        return copyOfRecord.build();
    }
}
