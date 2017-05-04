package org.campagnelab.dl.somatic.mappers.functional;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;
import java.util.Set;
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

    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across all samples and counts of these samples.
     * @param baseInformationOrBuilder
     * @param function
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forOneSampleGenotype(int sampleIndex,
                                                                                        int genotypeIndex,
                                                                                        BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                   Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> function) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();
        list.addAll(function.apply(baseInformationOrBuilder.getSamples(sampleIndex).getCounts(genotypeIndex)));
        return list;
    }

    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across all samples and counts of these samples, and both strands.
     * This CAN output a NumberWithFrequency list with repeat numbers. The subsequent mapper (like a DensityMapper) must handle this appropriately.
     * @param baseInformationOrBuilder
     * @param forwardFunction function to get NumberWithFrequency of forward strand
     * @param reverseFunction function to get NumberWithFrequency of reverse strand
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forOneSampleGenotypeBothStrands(int sampleIndex,
                                                                                        int genotypeIndex,
                                                                                        BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                        Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> forwardFunction,
                                                                                        Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> reverseFunction) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();
        BaseInformationRecords.CountInfo countInfo = baseInformationOrBuilder.getSamples(sampleIndex).getCounts(genotypeIndex);
        list.addAll(forwardFunction.apply(countInfo));
        list.addAll(reverseFunction.apply(countInfo));
        return list;
    }


    /**
     * Define a Function to reduce a record to a list of NumberWithFrequency found across N samples and counts of these samples.
     * @param baseInformationOrBuilder
     * @param function
     * @return
     */
    public static List<BaseInformationRecords.NumberWithFrequency> forNSampleCounts(Set<Integer> sampleIndices, BaseInformationRecords.BaseInformationOrBuilder baseInformationOrBuilder,
                                                                                      Function<BaseInformationRecords.CountInfo,List<BaseInformationRecords.NumberWithFrequency>> function) {
        List<BaseInformationRecords.NumberWithFrequency> list = new ObjectArrayList<>();

        for (int i = 0; i < baseInformationOrBuilder.getSamplesList().size(); i++) {
            if (!sampleIndices.contains(i)) {
                continue;
            }
            BaseInformationRecords.SampleInfo si = baseInformationOrBuilder.getSamples(i);
            for (BaseInformationRecords.CountInfo ci : si.getCountsList()){
                list.addAll(function.apply(ci));
            }
        }
        return list;
    }
}
