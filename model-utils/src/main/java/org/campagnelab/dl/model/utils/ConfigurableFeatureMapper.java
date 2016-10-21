package org.campagnelab.dl.model.utils;

import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;

/**
 * A configure method to customize mapping for properties of entire datasets.
 */
public interface ConfigurableFeatureMapper {
    /**
     * Configure the feature mapper for a specific set of sbi files. This method may access the properties of the reader
     * to retrieve statistics about the data being mapped to features (such as dataset wide normalization data).
     *
     * @param reader
     */
    void configure(SequenceBaseInformationReader reader);
}
