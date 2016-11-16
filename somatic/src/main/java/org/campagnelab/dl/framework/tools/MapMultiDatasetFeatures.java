package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Date;
import java.util.List;
import java.util.Properties;

/**
 * A tool to cache DL4J's multidatasets into in a .cf (cached features) file.
 *
 * @author Fabien Campagne
 */
public abstract class MapMultiDatasetFeatures<RecordType> extends AbstractTool<MapFeaturesArguments> {
    static private Logger LOG = LoggerFactory.getLogger(MapMultiDatasetFeatures.class);

    public void setNumRecordsWritten(int numRecordsWritten) {
        this.numRecordsWritten = numRecordsWritten;
    }

    protected abstract DomainDescriptor<RecordType> domainDescriptor();

    private int numRecordsWritten;

    @Override
    public MapFeaturesArguments createArguments() {
        return new MapMultiDatasetFeaturesArguments<RecordType>();
    }

    public MapMultiDatasetFeaturesArguments args() {
        return (MapMultiDatasetFeaturesArguments) arguments;
    }

    @Override
    public void execute() {
        //   assert args().adapter != null : "iterables must be provided in arguments.";
        DomainDescriptor<RecordType> domainDescriptor = domainDescriptor();
        MultiDataSetIteratorAdapter<RecordType> adapter = args().adapter;
        if (adapter == null) {

            try {
                adapter = new MultiDataSetIteratorAdapter<RecordType>(domainDescriptor.getRecordIterable(args().trainingSets,
                        (int) args().cacheN),
                        args().miniBatchSize,
                        domainDescriptor) {
                    @Override
                    public String getBasename() {
                        return buildBaseName(args().trainingSets);
                    }
                };
            } catch (IOException e) {
                throw new RuntimeException("Unable to load training set ", e);
            }
        }
       /* if ( != null) {
            inputs = Iterables.limit(Iterables.concat(args().iterables), (int) args().cacheN);
        } else {
            inputs = domainDescriptor.getRecordIterable(args().trainingSets, (int) args().cacheN);
        }
*/
           /* MultiDataSetIteratorAdapter<RecordType> adapter = new MultiDataSetIteratorAdapter<RecordType>(inputs, args().miniBatchSize, domainDescriptor) {
                @Override
                public String getBasename() {
                    return "multi-" + args().getTrainingSets().hashCode();
                }
            };

            LabelMapper labelMapper = new SimpleFeatureCalculator();
*/
        final String outputFilename = args().outputBasename + ".cf";
        try (FastBufferedOutputStream outputStream = new FastBufferedOutputStream(new FileOutputStream(outputFilename))) {
            ProgressLogger pg = new ProgressLogger(LOG);
            long numExamples = domainDescriptor.getNumRecords(args().getTrainingSets());
            pg.expectedUpdates = Math.min(numExamples, args().cacheN) / args().miniBatchSize;
            pg.displayLocalSpeed = true;
            pg.itemsName = "miniBatch";
            pg.start();
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            long numDatasets = 0;
            long writeAtMostN = args().writeAtMostN;
            numRecordsWritten = 0;

            while (adapter.hasNext()) {
                MultiDataSet mds = adapter.next();
                baos.reset();
                mds.save(baos);

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
                int numExamplesInDataset = mds.getFeatures()[0].size(0);
                numRecordsWritten += numExamplesInDataset;
                if (numRecordsWritten > args().cacheN) {
                    break;
                }
            }
            outputStream.close();
            pg.stop();

            long numRecords = domainDescriptor.getNumRecords(args().getTrainingSets());
            Properties cfpProperties = new Properties();
            cfpProperties.put("domainDescriptor", domainDescriptor().getClass().getCanonicalName());
            cfpProperties.put("multiDataSet", "true");
            cfpProperties.put("miniBatchSize", Integer.toString(args().miniBatchSize));
            if (args().domainDescriptor != null) {
                args().domainDescriptor.putProperties(cfpProperties);
            } else {
                cfpProperties.put("featureMapper", args().featureMapperClassname);
            }
            //  cfpProperties.put("labelMapper", labelMapper.getClass().getCanonicalName());
            cfpProperties.put("isTrio", Boolean.toString(args().isTrio));
            cfpProperties.put("numRecords", Long.toString(numRecordsWritten));
            cfpProperties.put("numDatasets", Long.toString(numDatasets));
            String[] inputNames = domainDescriptor.getComputationalGraph().getInputNames();
            for (String inputName : inputNames) {
                int dimIndex = 0;
                for (int dim : domainDescriptor().getNumInputs(inputName)) {
                    cfpProperties.put(inputName + ".numFeatures.dim" + Integer.toString(dimIndex), Integer.toString(dim));
                    dimIndex++;
                }
            }
            if (inputNames.length == 1) {
                // also write simpler numFeatures, for backward compatibility:
                cfpProperties.put("numFeatures", Integer.toString(domainDescriptor().getNumInputs(inputNames[0])[0]));
            }
            cfpProperties.put("stored", args().trainingSets.toString());
            cfpProperties.store(new FileWriter(new File(args().outputBasename + ".cfp")), new Date().toString());

        } catch (FileNotFoundException e) {
            LOG.error("Unable to create output file: " + outputFilename, e);
        } catch (IOException e) {
            LOG.error("Unable to write to output file: " + outputFilename, e);
        }

    }

    private String buildBaseName(List<String> trainingSets) {
        String cacheName;// only one input, use its name as cache name:
        if (trainingSets.size() == 1) {

            cacheName = FilenameUtils.getBaseName(trainingSets.get(0));
            ;
        } else {
            long hashcode = 8723872838723L;
            for (String name : trainingSets) {
                hashcode ^= FilenameUtils.getBaseName(name).hashCode();
            }
            cacheName = "multiset-" + Long.toString(hashcode);

        }
        return cacheName;
    }

    public int getNumRecordsWritten() {
        return numRecordsWritten;
    }

    public void setArguments(MapMultiDatasetFeaturesArguments<RecordType> arguments) {
        this.arguments = arguments;
    }
}
