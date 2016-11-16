package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.tools.MapMultiDatasetFeatures;
import org.campagnelab.dl.somatic.learning.domains.SomaticMutationDomainDescriptor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Properties;

/**
 * A tool to process an .sbi file and produce features in a .cf (cached features) file.
 *
 * @author Fabien Campagne
 */
public class MapFeaturesS extends MapMultiDatasetFeatures<BaseInformationRecords.BaseInformation> {
    static private Logger LOG = LoggerFactory.getLogger(MapFeaturesS.class);

    public static void main(String[] args) {

        MapFeaturesS tool = new MapFeaturesS();
        tool.parseArguments(args, "MapFeatures", tool.createArguments());
        tool.execute();
    }

    public void setNumRecordsWritten(int numRecordsWritten) {
        this.numRecordsWritten = numRecordsWritten;
    }

    @Override
    protected DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor() {
        Properties props = new Properties();
        props.put("input.featureMapper", args().featureMapperClassname);
        Properties sbiProperties = new Properties();
        // read statistics from the first training set:
        try (RecordReader reader = new RecordReader(args().getTrainingSets()[0])) {
            sbiProperties.putAll(reader.getProperties());
        } catch (IOException e) {
            throw new RuntimeException("Unable to obtain stats from training set " + args().getTrainingSets()[0],e);
        }
        return new SomaticMutationDomainDescriptor(props, sbiProperties);
    }

    private int numRecordsWritten;

    @Override
    public MapFeaturesArguments createArguments() {
        return new MapFeaturesArguments();
    }

    public void setArguments(MapFeaturesArguments arguments) {
        this.arguments = arguments;
    }
}
