package org.campagnelab.dl.somatic.mappers;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by fac2003 on 12/21/16.
 */
public class MappingFunctions {
    public static String recordTo(final int contextLength, BaseInformationRecords.BaseInformationOrBuilder record, int countIndex) {
        final List<BaseInformationRecords.CountInfo> counts = record.getSamples(record.getSamplesCount() - 1).getCountsList();
        ArrayList<BaseInformationRecords.CountInfo> sorted = new ArrayList<>();
        sorted.addAll(counts);
        Collections.sort(sorted, (o1, o2) ->
                (o2.getGenotypeCountForwardStrand() + o2.getGenotypeCountReverseStrand()) - (o1.getGenotypeCountForwardStrand() + o1.getGenotypeCountReverseStrand())
        );
        String to = countIndex < sorted.size() ? sorted.get(countIndex).getToSequence() : "N";
        StringBuffer context = new StringBuffer(to);

        while (context.length() < contextLength) {
            context.append("N");
        }
        context.setLength(contextLength);
        return context.toString();
    }
}
