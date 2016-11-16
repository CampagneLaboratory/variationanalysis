package org.campagnelab.dl.somatic.intermediaries;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Mutate the second sample in a pair of sample to simulate somatic mutation and errors.
 * Created by fac2003 on 7/19/16.
 */
public interface SimulationStrategy {
    BaseInformationRecords.BaseInformation mutate(boolean makeSomatic,
                                                  BaseInformationRecords.BaseInformation record,
                                                  BaseInformationRecords.SampleInfo germlineSample,
                                                  BaseInformationRecords.SampleInfo otherSample, SimulationCharacteristics sim);

    void setSeed(long seed);
}
