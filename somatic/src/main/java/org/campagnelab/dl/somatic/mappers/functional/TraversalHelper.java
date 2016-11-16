package org.campagnelab.dl.somatic.mappers.functional;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;
import java.util.function.Function;

/**
 * Created by fac2003 on 10/21/16.
 */
public class TraversalHelper {

    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across all samples and counts of these samples.
     * @param baseInformationOrBuilder
     * @param function
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forAllSampleCounts(BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> function) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();

        baseInformationOrBuilder.getSamplesList().forEach(
                sampleInfo -> {
                    sampleInfo.getCountsList().forEach(
                            countInfo -> list.addAll(function.apply(countInfo))
                    );
                }
        );
        return list;
    }

    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across all samples and counts of these samples.
     * @param baseInformationOrBuilder
     * @param function
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forSampleCounts(int sampleIndex, BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> function) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();

        baseInformationOrBuilder.getSamples(sampleIndex).getCountsList().forEach(
                            countInfo -> list.addAll(function.apply(countInfo))
                    );
        return list;
    }
}
