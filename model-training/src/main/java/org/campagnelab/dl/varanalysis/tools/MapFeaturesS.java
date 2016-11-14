package org.campagnelab.dl.varanalysis.tools;

import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.learning.TrainSomaticModel;
import org.campagnelab.dl.varanalysis.learning.domains.DomainDescriptor;
import org.campagnelab.dl.varanalysis.learning.domains.SomaticMutationDomainDescriptor;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Date;
import java.util.List;
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
