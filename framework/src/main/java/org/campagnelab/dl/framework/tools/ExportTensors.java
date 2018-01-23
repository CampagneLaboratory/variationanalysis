package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapterMultipleSamples;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.datasets.iterator.AsyncMultiDataSetIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

/**
 * A tool to export mapped features into tensor files. Useful to write files that can be loaded as numpy tensors in python.
 *
 * @author Fabien Campagne
 */
public abstract class ExportTensors<RecordType> extends AbstractTool<ExportTensorArguments> {
    static private Logger LOG = LoggerFactory.getLogger(ExportTensors.class);

    protected abstract DomainDescriptor<RecordType> domainDescriptor(String featureMapperClassName,
                                                                     List<String> trainingSets);

    protected abstract void decorateProperties(Properties decorateProperties);

    private int numRecordsWritten;

    @Override
    public ExportTensorArguments createArguments() {
        return new ExportTensorArguments();
    }

    public ExportTensorArguments args() {
        return (ExportTensorArguments) arguments;
    }

    @Override
    public void execute() {
        if ((args().sampleNames.size() != args().sampleTypes.size())
                || (args().sampleNames.size() != args().sampleIds.size())) {
            throw new RuntimeException("--sample-name,  --sample-type and --sample-index must be provided the same " +
                    "number of times on the command line.");
        }
        DomainDescriptor<RecordType> domainDescriptor = domainDescriptor(args().featureMapperClassname,
                args().trainingSets);
        MultiDataSetIteratorAdapterMultipleSamples<RecordType> adapter;
        IntArrayList sampleIndices = new IntArrayList();
        sampleIndices.addAll(args().sampleIds);
        long numRecordsLong = domainDescriptor.getNumRecords(args().getTrainingSets());
        if (numRecordsLong > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Too many records to write in training sets");
        }
        int numRecords = Integer.min(args().exportN, (int) numRecordsLong);


        try {
            Properties properties = getReaderProperties(args().trainingSets.get(0));
            decorateProperties(properties);
            adapter = new MultiDataSetIteratorAdapterMultipleSamples<RecordType>(domainDescriptor.getRecordIterable(args().trainingSets,
                    (int) args().exportN),
                    args().miniBatchSize,
                    domainDescriptor,
                    false,
                    null,
                    sampleIndices.toIntArray(),
                    properties) {
                @Override
                public String getBasename() {
                    return buildBaseName(args().trainingSets);
                }
            };
        } catch (IOException e) {
            throw new RuntimeException("Unable to load training set ", e);
        }

        Iterator<List<MultiDataSet>> iterator = adapter;
        ProgressLogger pg = new ProgressLogger(LOG);
        long numExamples = domainDescriptor.getNumRecords(args().getTrainingSets());
        pg.expectedUpdates = Math.min(numExamples, args().exportN) / args().miniBatchSize;
        pg.displayLocalSpeed = true;
        pg.itemsName = "miniBatch";
        pg.start();
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        long numDatasets = 0;
        numRecordsWritten = 0;
        String[] inputNames = domainDescriptor.getComputationalGraph().getInputNames();
        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();

        System.out.println("All inputs: " + ObjectArrayList.wrap(inputNames));
        System.out.println("All outputs: " + ObjectArrayList.wrap(outputNames));
        final Set<String> inputNamesToExport = args().inputNamesToExport;
        final Set<String> outputNamesToExport = args().outputNamesToExport;

        int[] inputIndicesSelected = assignSelectedIndices(inputNames, inputNamesToExport);
        int[] outputIndicesSelected = assignSelectedIndices(outputNames, outputNamesToExport);


        int miniBatchSize = args().miniBatchSize;
        try (final VectorWriter vectorWriter = new VectorWriterText(args().outputBasename)) {

            // TODO: Set from args or properties
            vectorWriter.setSpecVersionNumber(0, 1);
            vectorWriter.setNumRecords(numRecords);

            int currStartExampleIndex = 0;
            while (iterator.hasNext()) {
                List<MultiDataSet> mdsList = iterator.next();
                vectorWriter.appendMdsList(mdsList, inputIndicesSelected, outputIndicesSelected, inputNames,
                        outputNames, currStartExampleIndex);
                currStartExampleIndex += miniBatchSize;
                pg.update();
            }

            vectorWriter.addSampleInfo("testSampleType", "testSampleName");
        } catch (IOException e) {
            throw new RuntimeException("Unable to write vector file. ", e);
        } finally {
            pg.stop();
        }

        Properties cfpProperties = new Properties();
        cfpProperties.put("domainDescriptor", domainDescriptor.getClass().getCanonicalName());
        cfpProperties.put("multiDataSet", "true");
        cfpProperties.put("miniBatchSize", Integer.toString(args().miniBatchSize));
        if (domainDescriptor != null) {
            domainDescriptor.putProperties(cfpProperties);
        } else {
            cfpProperties.put("featureMapper", args().featureMapperClassname);
        }
        cfpProperties.put("numRecords", Long.toString(numRecordsWritten));
        cfpProperties.put("numDatasets", Long.toString(numDatasets));
        for (String inputName : inputNames) {
            int dimIndex = 0;
            for (int dim : domainDescriptor.getNumInputs(inputName)) {
                cfpProperties.put(inputName + ".numFeatures.dim" + Integer.toString(dimIndex), Integer.toString(dim));
                dimIndex++;
            }
        }
        if (inputNames.length == 1) {
            // also write simpler numFeatures, for backward compatibility:
            cfpProperties.put("numFeatures", Integer.toString(domainDescriptor.getNumInputs(inputNames[0])[0]));
        }
        cfpProperties.put("stored", args().trainingSets.toString());
        try (final FileWriter writer = new FileWriter(new File(args().outputBasename + ".cfp"))) {
            cfpProperties.store(writer, new Date().toString());
        } catch (IOException e) {
            throw new RuntimeException("Unable to write .vecp file.", e);
        }

    }

    private int[] assignSelectedIndices(String[] inputNames, Set<String> inputNamesToExport) {
        if (inputNamesToExport == null) {
            inputNamesToExport = new ObjectArraySet<>();
            inputNamesToExport.addAll(ObjectArrayList.wrap(inputNames));
        }
        int[] inputIndicesSelected = new int[inputNames.length];
        if (!inputNamesToExport.isEmpty()) {
            inputIndicesSelected = new int[inputNamesToExport.size()];
            Arrays.fill(inputIndicesSelected, -1);
            int i = 0;
            for (int j = 0; j < inputNames.length; j++) {
                if (inputNamesToExport.contains(inputNames[j])) {
                    inputIndicesSelected[i++] = j;
                }
            }
            for (int index : inputIndicesSelected) {
                if (index == -1) {
                    System.err.printf("An argument to --export-inputs/outputs does not match actual inputs/outputs.");
                    System.exit(1);
                }
            }

        }
        return inputIndicesSelected;
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


    // Taken from org.campagnelab.dl.somatic.learning.SomaticTrainer
    private Properties getReaderProperties(String trainingSet) throws IOException {
        try (SequenceBaseInformationReader reader = new SequenceBaseInformationReader(trainingSet)) {
            final Properties properties = reader.getProperties();
            reader.close();
            return properties;
        }
    }
}
