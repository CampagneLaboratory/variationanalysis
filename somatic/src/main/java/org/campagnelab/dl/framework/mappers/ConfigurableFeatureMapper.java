package org.campagnelab.dl.framework.mappers;


import java.util.Properties;

/**
 * A configure method to customize mapping for properties of entire datasets.
 */
public interface ConfigurableFeatureMapper {
    /**
     * Configure the feature mapper for a specific set of sbi files. This method may access the properties of the reader
     * to retrieve statistics about the data being mapped to features (such as dataset wide normalization data).
     *
     * @param readerProperties properties from a goby sbi reader.
     */
    void configure(Properties readerProperties);
}
