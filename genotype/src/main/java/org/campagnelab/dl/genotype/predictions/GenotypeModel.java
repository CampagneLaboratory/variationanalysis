package org.campagnelab.dl.genotype.predictions;

import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.DomainDescriptorLoader;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.algorithmic.dsv.DiscoverVariantPositionData;
import org.campagnelab.goby.algorithmic.dsv.SampleCountInfo;
import org.campagnelab.goby.predictions.ProtoHelper;
import org.campagnelab.goby.reads.RandomAccessSequenceInterface;
import org.deeplearning4j.nn.api.Model;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Constructor;
import java.util.Properties;

/**
 * Encapsulates a trained genotype model.
 * Created by fac2003 on 12/18/16.
 */
public class GenotypeModel {
    static private Logger LOG = LoggerFactory.getLogger(GenotypeModel.class);
    private final String modelPath;
    private final Properties modelProperties;
    private DomainDescriptor domainDescriptor;
    private GenotypeProtoPredictor protoPredictor;

    /**
     * Create a GenotypeModel with model path and prefix/label.
     */
    public GenotypeModel(String modelPath, String prefix) throws IOException {
        this.modelPath = modelPath;
        //get MAPPER
        FeatureMapper featureMapper = null;
        Properties prop = new Properties();
        InputStream input = null;
        String mapperName = null;
        try {
            final String modelPropertiesFilename = modelPath + "/config.properties";
            if (!new File(modelPropertiesFilename).exists()) {
                LOG.warn("model property file does not exist: " + modelPropertiesFilename);
            }
            input = new FileInputStream(modelPropertiesFilename);
            // load a properties file
            prop.load(input);
            // get the property value and print it out
            mapperName = prop.getProperty("mapper");
            if (mapperName == null) {
                LOG.warn("property mapper in model config.properties file is not defined.");
            }
            ClassLoader classLoader = this.getClass().getClassLoader();
            // Load the target class using its binary name
            Class loadedMyClass = classLoader.loadClass(mapperName);
            System.out.println("Loaded class name: " + loadedMyClass.getName());
            // Create a new instance from the loaded class
            Constructor constructor = loadedMyClass.getConstructor();
            featureMapper = (FeatureMapper) constructor.newInstance();
            if (featureMapper instanceof ConfigurableFeatureMapper) {
                ConfigurableFeatureMapper confMapper = (ConfigurableFeatureMapper) featureMapper;
                System.out.println("Configuring feature mapper with model properties at " + modelPropertiesFilename);
                confMapper.configure(prop);
            }
        } catch (Exception e) {
            throw new RuntimeException("Unable to create feature mapper " + mapperName, e);
        } finally {
            IOUtils.closeQuietly(input);
        }


        ModelLoader modelLoader = new ModelLoader(modelPath);
        Model model = modelLoader.loadModel(prefix);
        modelProperties = modelLoader.getModelProperties();
        domainDescriptor = DomainDescriptorLoader.load(modelPath);
        this.protoPredictor = new GenotypeProtoPredictor(domainDescriptor, model, featureMapper);
    }

    public GenotypePrediction predictGenotype(RandomAccessSequenceInterface genome, String referenceID,
                                                  SampleCountInfo sampleCounts[],
                                                  int referenceIndex, int position,
                                                  DiscoverVariantPositionData list,
                                                  int[] readerIdxs) {
        Integer[] sampleToReaderIdxs;
        // genotype models work with a single sample:
        sampleToReaderIdxs = new Integer[]{readerIdxs[0]};

        //in the past, predictions on 0 reads have been bypassed and given prediction value 0. leaving this out for now.
        BaseInformationRecords.BaseInformation proto = ProtoHelper.toProto(genome, referenceID, sampleCounts, referenceIndex, position, list, sampleToReaderIdxs);
        return protoPredictor.predictGenotype(proto);
    }


    public Properties getProperties() {
        return modelProperties;
    }
}
