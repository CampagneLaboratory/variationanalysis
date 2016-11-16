package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.tools.MapFeaturesArguments;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.nd4j.linalg.dataset.DataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Date;
import java.util.Properties;

/**
 * A tool to process an .sbi file and produce features in a .cf (cached features) file.
 *
 * @author Fabien Campagne
 */
public class MapFeatures extends AbstractTool<SomaticMapFeaturesArguments> {
    static private Logger LOG = LoggerFactory.getLogger(MapFeatures.class);

    public void setNumRecordsWritten(int numRecordsWritten) {
        this.numRecordsWritten = numRecordsWritten;
    }

    public void setArguments(SomaticMapFeaturesArguments arguments) {
        this.arguments = arguments;
    }

    private int numRecordsWritten;

    public static void main(String[] args) {

        MapFeatures tool = new MapFeatures();
        tool.parseArguments(args, "MapFeatures", tool.createArguments());
        tool.execute();
    }

    @Override
    public SomaticMapFeaturesArguments createArguments() {
        return new SomaticMapFeaturesArguments();
    }

    @Override
    public void execute() {
        BaseInformationConcatIterator inputs = null;
        LabelMapper labelMapper = new SimpleFeatureCalculator();
        try {
            String[] trainingDataset = args().getTrainingSets();
            FeatureMapper featureMapper =
                    TrainSomaticModel.configureFeatureMapper(args().featureMapperClassname, args().isTrio, trainingDataset);
            if (args().iterators != null) {
                inputs = new BaseInformationConcatIterator(args().iterators, args().miniBatchSize, featureMapper, labelMapper);
            } else {
                inputs = new BaseInformationConcatIterator(args().miniBatchSize, featureMapper, labelMapper, args().getTrainingSets());
            }
        } catch (IOException e) {
            LOG.error("Unable to load input files.");
            System.exit(1);
        }
        final String outputFilename = args().outputBasename + ".cf";
        try (FastBufferedOutputStream outputStream = new FastBufferedOutputStream(new FileOutputStream(outputFilename))) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.expectedUpdates = inputs.numExamples() / args().miniBatchSize;
            pg.displayLocalSpeed = true;
            pg.itemsName = "miniBatch";
            pg.start();
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            long numDatasets = 0;
            long writeAtMostN = args().writeAtMostN;
            numRecordsWritten = 0;

            while (inputs.hasNext()) {
                DataSet ds = inputs.next();
                baos.reset();
                ds.save(baos);

                final byte[] bytes = baos.toByteArray();

                // write the length of the array first, most significant bytes first:
                outputStream.write((bytes.length >> 8 * 3) & 0xFF);
                outputStream.write((bytes.length >> 8 * 2) & 0xFF);
                outputStream.write((bytes.length >> 8) & 0xFF);
                outputStream.write(bytes.length & 0xFF);
                outputStream.write(bytes);
                pg.lightUpdate();
                if (numRecordsWritten > writeAtMostN) {
                    break;
                }
                numDatasets += 1;
                numRecordsWritten += ds.numExamples();
                if (numRecordsWritten > args().cacheN) {
                    break;
                }
            }
            outputStream.close();
            pg.stop();

            Properties cfpProperties = new Properties();
            long numRecords = 0;
            for (String inputFile : args().getTrainingSets()) {
                SequenceBaseInformationReader reader = new SequenceBaseInformationReader(inputFile);
                final Properties sbiProperties = reader.getProperties();
                cfpProperties.putAll(sbiProperties);
                numRecords += Long.parseLong(sbiProperties.get("numRecords").toString());
            }
            cfpProperties.put("multiDataSet", "false");
            cfpProperties.put("miniBatchSize", Integer.toString(args().miniBatchSize));
            cfpProperties.put("featureMapper", args().featureMapperClassname);
            cfpProperties.put("labelMapper", labelMapper.getClass().getCanonicalName());
            cfpProperties.put("isTrio", Boolean.toString(args().isTrio));
            cfpProperties.put("numRecords", Long.toString(numRecordsWritten));
            cfpProperties.put("numDatasets", Long.toString(numDatasets));
            cfpProperties.put("numFeatures", Integer.toString(inputs.inputColumns()));
            cfpProperties.put("stored", args().trainingSets.toString());
            cfpProperties.store(new FileWriter(new File(args().outputBasename + ".cfp")), new Date().toString());

        } catch (FileNotFoundException e) {
            LOG.error("Unable to create output file: " + outputFilename);
        } catch (IOException e) {
            LOG.error("Unable to write to output file: " + outputFilename);
        }
    }


    public int getNumRecordsWritten() {
        return numRecordsWritten;
    }

}
