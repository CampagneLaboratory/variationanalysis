package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapterMultipleSamples;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.nd4j.linalg.dataset.api.MultiDataSet;
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
                                                                     List<String> trainingSets,
                                                                     int genomicContextLength,
                                                                     float labelSmoothingEpsilon,
                                                                     int ploidy, int extraGenotypes);

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
                args().trainingSets, args().genomicContextLength, args().labelSmoothingEpsilon,
                args().ploidy, args().extraGenotypes);

        MultiDataSetIteratorAdapterMultipleSamples<RecordType> adapter;
        IntArrayList sampleIndices = new IntArrayList();
        sampleIndices.addAll(args().sampleIds);
        long numRecordsLong = Long.min(args().exportN, domainDescriptor.getNumRecords(args().getTrainingSets()));
        if (numRecordsLong > Integer.MAX_VALUE) {
            LOG.warn(String.format("Number of records to be written %d is greater than max possible value of %d and" +
                    " will be truncated", numRecordsLong, Integer.MAX_VALUE));
        }
        int numRecords = (int) numRecordsLong;


        try {
            Properties properties;
            if (args().sbiList == null) {
                properties = getReaderProperties(args().trainingSets.get(0));
            } else {
                properties = new Properties();
                String mergedSbipBaseName = FilenameUtils.removeExtension(args().sbiList) + "_merged";
                try {
                    BufferedReader sbiListReader = new BufferedReader(new FileReader(args().sbiList));
                    String currSbiPath;
                    ObjectList<Properties> sbiProps = new ObjectArrayList<>();
                    while ((currSbiPath = sbiListReader.readLine()) != null) {
                        sbiProps.add(getReaderProperties(currSbiPath));
                    }
                    SequenceBaseInformationWriter.writeProperties(mergedSbipBaseName, sbiProps);
                    Reader propertiesReader = new FileReader(mergedSbipBaseName + ".sbip");
                    properties.load(propertiesReader);
                    IOUtils.closeQuietly(sbiListReader);
                    IOUtils.closeQuietly(propertiesReader);
                } catch (IOException e) {
                    throw new RuntimeException("Unable to load combined sbips", e);
                }
            }
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
        try (final VectorWriter vectorWriter = createVectorWriter()) {
            vectorWriter.setNumRecords(numRecords);
            vectorWriter.setDomainDescriptor(domainDescriptor.getClass().getCanonicalName());
            vectorWriter.setFeatureMapper(args().featureMapperClassname);
            String[] trainingSets = new String[args().trainingSets.size()];
            trainingSets = args().trainingSets.toArray(trainingSets);
            vectorWriter.setInputFiles(trainingSets);
            long currStartExampleIndex = 0;
            while (iterator.hasNext()) {
                List<MultiDataSet> mdsList = iterator.next();
                vectorWriter.appendMdsList(mdsList, inputIndicesSelected, outputIndicesSelected, inputNames,
                        outputNames, currStartExampleIndex);
                currStartExampleIndex += miniBatchSize;
                pg.update();
            }
            for (int sampleIndex : args().sampleIds) {
                vectorWriter.addSampleInfo(args().sampleTypes.get(sampleIndex), args().sampleNames.get(sampleIndex));
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to write vector file. ", e);
        } finally {
            pg.stop();
        }
        try {
            String path=FilenameUtils.getFullPath(args().outputBasename);
            if (path.equals("/")) {
                path="/.";
            }
            String inputFiles = Arrays.toString(args().getTrainingSets());
            Properties propertiesToExport = getReaderProperties(args().trainingSets.get(0));
            decorateProperties(propertiesToExport);
            domainDescriptor.writeProperties(path, propertiesToExport);
            Properties configPropertiesToExport = new Properties();
            String configPropertiesPath = ModelLoader.getModelPath(path) + "/config.properties";
            configPropertiesToExport.putAll(propertiesToExport);
            configPropertiesToExport.put("domainDescriptor", domainDescriptor.getClass().getCanonicalName());
            configPropertiesToExport.store(new FileWriter(configPropertiesPath),
                    "Config properties created from export tensors on files " + inputFiles);
        } catch (IOException e) {
            System.err.println("Unable to write domain descriptor properties");
        }
    }

    private VectorWriter createVectorWriter() throws IOException {
        VectorWriter vectorWriter;
        if (args().outputBasename==null && args().trainingSets.size()==1) {
            // construct the output basename with a path:
            final String filename = args().trainingSets.iterator().next();
            String inputBasename=FilenameUtils.getPath(filename)+FilenameUtils.getBaseName(filename);
            if (filename.startsWith("/")) {
                inputBasename="/"+inputBasename;
            }
            if (FilenameUtils.getPathNoEndSeparator(filename).equals("")){
                inputBasename="./"+inputBasename;
            }
            args().outputBasename=inputBasename;
        }
        if (args().outputBasename==null){
            System.err.println("Please specify -o when using multiple inputs.");
            System.exit(1);
        }
        switch (args().vecFileType) {
            case "text":
            case "gzipped+text":
                vectorWriter = new VectorWriterText(args().outputBasename);
                break;
            case "binary":
                vectorWriter = new VectorWriterBinary(args().outputBasename);
                break;
            default:
                throw new IllegalArgumentException(
                        String.format("Unknown vector file type %s; choices are text and binary.",
                                args().vecFileType)
                );
        }
        return vectorWriter;
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
    public abstract Properties getReaderProperties(String trainingSet) throws IOException;



}
