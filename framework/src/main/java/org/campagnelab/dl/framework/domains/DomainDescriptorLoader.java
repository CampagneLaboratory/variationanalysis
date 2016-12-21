package org.campagnelab.dl.framework.domains;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.util.Properties;

/**
 * Helper to load the domain descriptor associated with a set of models.
 * Created by fac2003 on 11/12/16.
 */
public class DomainDescriptorLoader {
    static private Logger LOG = LoggerFactory.getLogger(DomainDescriptorLoader.class);

    /**
     * Load the model descriptors associated with models in the model directory. A runtime
     * exception is thrown is any errors occurs at this point.
     *
     * @param modelPath model director.
     * @return A model descriptor.
     */
    public static DomainDescriptor load(String modelPath) {
        try {
            Properties modelProperties = new Properties();
            String configFilename = modelPath + "/config.properties";
            FileInputStream input = new FileInputStream(configFilename);
            // load a properties file

            modelProperties.load(input);
            // get the property value and print it out
            String domainDescriptorClassname = modelProperties.getProperty("domainDescriptor");
           /* if (domainDescriptorClassname == null) {
                LOG.warn("domainDescriptor property was not found in " + configFilename);
                LOG.warn("Using backward compatible descriptor instead: SomaticMutationDomainDescriptor");
                return new SomaticMutationDomainDescriptor(modelPath);
            }*/
            // we call the constructor that accepts a model path as single argument:
            Class<DomainDescriptor> clazz =
                    (Class<DomainDescriptor>) Class.forName(domainDescriptorClassname);
            final DomainDescriptor domainDescriptor = clazz.getConstructor(String.class).newInstance(modelPath);
            domainDescriptor.configure(modelProperties);
            return domainDescriptor;

        } catch (Exception e) {
            throw new RuntimeException("Unable to load domain descriptor with model path: " +modelPath,e);
        }
    }
}
