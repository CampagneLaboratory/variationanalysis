package org.campagnelab.dl.framework.domains;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
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

    /**
     * Load the model descriptors associated with a properties file and a specified class. A runtime
     * exception is thrown is any errors occurs at this point.
     *
     * @param domainPropertiesPath path to domain.properties file
     * @param domainClass fully qualified classname for domain descriptor implementation
     * @return A model descriptor.
     */
    public static DomainDescriptor loadFromProperties(String domainPropertiesPath,
                                                      String domainClass) {
        DomainDescriptor domainDescriptor;
        try (FileInputStream input = new FileInputStream(domainPropertiesPath)) {
            Properties domainProperties = new Properties();
            domainProperties.load(input);
            Class<DomainDescriptor> clazz = (Class<DomainDescriptor>) Class.forName(domainClass);
            domainDescriptor = clazz.getConstructor(Properties.class, Properties.class).newInstance(domainProperties,
                    new Properties());
        } catch (IOException e) {
            throw new RuntimeException("Couldn't load domain descriptor from " + domainPropertiesPath, e);
        } catch (ClassNotFoundException | NoSuchMethodException | InstantiationException | IllegalAccessException
                | InvocationTargetException e) {
            throw new RuntimeException("Couldn't instantiate " + domainClass, e);
        }
        return domainDescriptor;
    }
}
