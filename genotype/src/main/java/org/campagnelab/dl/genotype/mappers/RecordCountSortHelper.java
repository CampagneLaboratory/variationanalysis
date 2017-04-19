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
        return sort(0, record);
    }

    public BaseInformationRecords.BaseInformationOrBuilder sort(int sampleIndex, BaseInformationRecords.BaseInformationOrBuilder record) {
        int originalGenotypeIndex = 0;
        final List<BaseInformationRecords.CountInfo> countsList = record.getSamples(sampleIndex).getCountsList();
        List<BaseInformationRecords.CountInfo> counts = new ArrayList<>();
        for (BaseInformationRecords.CountInfo count : countsList) {
            counts.add(BaseInformationRecords.CountInfo.newBuilder()
                    .mergeFrom(count)
                    .setGobyGenotypeIndex(originalGenotypeIndex++)
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
        final BaseInformationRecords.SampleInfo.Builder builder = record.getSamples(sampleIndex).toBuilder();
        builder.clearCounts();
        builder.addAllCounts(counts);
        copyOfRecord.addAllSamples(record.getSamplesList());
       // overwrite the sample we have sorted:
        copyOfRecord.setSamples(sampleIndex, builder);
        return copyOfRecord.build();
    }

    public BaseInformationRecords.BaseInformation sort(BaseInformationRecords.BaseInformation record) {
        return sort(0, record);
    }

    public BaseInformationRecords.BaseInformation sort(int sampleIndex, BaseInformationRecords.BaseInformation record) {
        int originalGenotypeIndex = 0;
        final List<BaseInformationRecords.CountInfo> countsList = record.getSamples(sampleIndex).getCountsList();
        List<BaseInformationRecords.CountInfo> counts = new ArrayList<>();
        for (BaseInformationRecords.CountInfo count : countsList) {
            counts.add(BaseInformationRecords.CountInfo.newBuilder()
                    .mergeFrom(count)
                    .setGobyGenotypeIndex(originalGenotypeIndex++)
                    .build());
        }

        Collections.sort(counts, (a, b) ->
                (b.getGenotypeCountForwardStrand() + b.getGenotypeCountReverseStrand()) -
                        (a.getGenotypeCountForwardStrand() + a.getGenotypeCountReverseStrand()));


        final BaseInformationRecords.SampleInfo.Builder builder = record.getSamples(sampleIndex).toBuilder();
        builder.clearCounts();
        builder.addAllCounts(counts);
        BaseInformationRecords.BaseInformation.Builder copyOfRecord = record.toBuilder();
        copyOfRecord.setSamples(sampleIndex, builder);
        return copyOfRecord.build();
    }
}
