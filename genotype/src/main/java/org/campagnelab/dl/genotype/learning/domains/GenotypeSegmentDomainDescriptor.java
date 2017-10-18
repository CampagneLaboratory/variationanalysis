package org.campagnelab.dl.genotype.learning.domains;

import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.ConfigurableLabelMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousInterpreter;
import org.campagnelab.dl.genotype.mappers.NumDistinctAllelesLabelMapper;
import org.campagnelab.dl.genotype.mappers.SingleBaseFeatureMapperV1;
import org.campagnelab.dl.genotype.mappers.SingleBaseLabelMapperV1;
import org.campagnelab.dl.genotype.tools.SegmentTrainingArguments;
import org.campagnelab.dl.genotype.storage.SegmentReader;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationReader;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.lossfunctions.ILossFunction;
import org.nd4j.linalg.lossfunctions.impl.LossMCXENT;

import java.io.IOException;
import java.util.Map;
import java.util.Properties;
import java.util.function.Function;

/**
 * A domain descriptor for the task of predicting genotypes for segments/blocks of the genome.
 */
public class GenotypeSegmentDomainDescriptor extends DomainDescriptor<SegmentInformationRecords.SegmentInformation> {
    private final SegmentTrainingArguments arguments;
    private int ploidy;

    public GenotypeSegmentDomainDescriptor(SegmentTrainingArguments arguments) {
        this.arguments = arguments;
        initializeArchitecture(arguments.architectureClassname);
        this.ploidy = arguments.ploidy;
    }

    @Override
    public void putProperties(Properties props) {
        super.putProperties(props);
        props.setProperty(NumDistinctAllelesLabelMapper.PLOIDY_PROPERTY, Integer.toString(ploidy));
        decorateProperties(props);
    }

    /**
     * Record arguments to the properties, that need to be provided to feature/label mappers.
     *
     * @param modelProperties
     */
    void decorateProperties(Properties modelProperties) {
        // transfer arguments to properties when we know we started from command line arguments, otherwise the
        // arguments are already in properties:
        if (args().parsedFromCommandLine) {

            modelProperties.setProperty("stats.genomicContextSize.min", "1");
            modelProperties.setProperty("stats.genomicContextSize.max", "1");
            modelProperties.setProperty("indelSequenceLength", "1");


            ploidy = args().ploidy;
            modelProperties.setProperty("genotypes.ploidy", Integer.toString(args().ploidy));
            modelProperties.setProperty("genotypes.ploidy", Integer.toString(args().ploidy));
        }
    }

    private SegmentTrainingArguments args() {
        return arguments;
    }

    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        ConfigurableFeatureMapper mapper;
        if (cachedFeatureMappers.containsKey(inputName)) {
            return cachedFeatureMappers.get(inputName);
        } else {
            switch (inputName) {
                case "input":
                    mapper = new SingleBaseFeatureMapperV1(0);
                    break;
                default:
                    throw new RuntimeException("Unsupported input name: " + inputName);
            }
            if (mapper != null) {
                try {

                    Properties properties = getReaderProperties(args().trainingSets.get(0));
                    decorateProperties(properties);
                    mapper.configure(properties);
                    cachedFeatureMappers.put(inputName, (FeatureMapper) mapper);
                    return (FeatureMapper) mapper;
                } catch (IOException e) {
                    throw new RuntimeException("IO exception, perhaps .sbip file not found?", e);
                }

            } else {
                return null;

            }

        }
    }
    public static Properties getReaderProperties(String trainingSet) throws IOException {
        try (SequenceSegmentInformationReader reader = new SequenceSegmentInformationReader(trainingSet)) {
            final Properties properties = reader.getProperties();
            reader.close();
            return properties;
        }
    }

    Map<String, FeatureMapper> cachedFeatureMappers = new Object2ObjectOpenHashMap<>();
    Map<String, LabelMapper> cachedLabelMappers = new Object2ObjectOpenHashMap<>();


    @Override
    public LabelMapper getLabelMapper(String outputName) {
        if (cachedLabelMappers.containsKey(outputName)) {
            return cachedLabelMappers.get(outputName);
        } else {
            switch (outputName) {
                case "genotype":
                    try {
                        ConfigurableLabelMapper mapper = new SingleBaseLabelMapperV1(0);
                        final Properties readerProperties = getReaderProperties(args().trainingSets.get(0));
                        decorateProperties(readerProperties);
                        mapper.configure(readerProperties);
                        cachedLabelMappers.put(outputName,(LabelMapper) mapper);
                        return (LabelMapper) mapper;
                    } catch(IOException e) {
                        throw new InternalError("Unable to load properties and initialize label mapper.",e);
                    }
                default:
                    throw new RuntimeException("Unsupported output name: " + outputName);
            }
        }
    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {
        // TODO: change this to the correct interpreter for this domain
        return new HomozygousInterpreter(false);
    }

    @Override
    public Function<String, ? extends Iterable<SegmentInformationRecords.SegmentInformation>> getRecordIterable() {
        return inputFilename -> {
            try {
                return new SegmentReader(inputFilename);
            } catch (IOException e) {
                throw new RuntimeException("Unable to read records from " + inputFilename, e);
            }
        };
    }

    private ComputationGraphAssembler assembler = null;

    @Override
    public synchronized ComputationGraphAssembler getComputationalGraph() {
        if (assembler == null) {
            assembler = new GenotypeSegmentsLSTM();
            assembler.setArguments(arguments);
        }
        return assembler;
    }

    @Override
    public int[] getNumInputs(String inputName) {
        return getFeatureMapper(inputName).dimensions().dimensions;
    }

    @Override
    public int[] getNumOutputs(String outputName) {
        return getLabelMapper(outputName).dimensions().dimensions;
    }

    @Override
    public int[] getNumMaskInputs(String inputName) {
        // drop the last dimension:
        int[] dimensions=new int[1];
        dimensions[0]=getNumInputs(inputName)[1];
        return dimensions;
    }

    @Override
    public int[] getNumMaskOutputs(String outputName) {
        // drop the last dimension:
        int[] dimensions=new int[1];
        dimensions[0]=getNumOutputs(outputName)[1];
        return dimensions;
    }

    @Override
    public int getNumHiddenNodes(String componentName) {
        return args().numHiddenNodes;
    }

    @Override
    public ILossFunction getOutputLoss(String outputName) {
        return new LossMCXENT();
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        long count = 0;
        for (String recordFile : recordFiles) {

            try (SegmentReader reader = new SegmentReader(recordFile)) {
                count += reader.getTotalRecords();
            } catch (IOException e) {
                return 0;
            }
        }
        return count;
    }
}
