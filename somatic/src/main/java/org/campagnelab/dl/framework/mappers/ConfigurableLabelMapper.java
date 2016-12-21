package org.campagnelab.dl.framework.mappers;

import java.util.Properties;

/**
 * A configure method to customize mapping for dataset labels
 */
public interface ConfigurableLabelMapper {
    /**
     * Configure the label mapper for a specific dataset.
     *
     * @param readerProperties properties containing relevant information
     */
    void configure(Properties readerProperties);
}
