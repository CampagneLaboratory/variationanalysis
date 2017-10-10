package org.campagnelab.dl.genotype.learning.domains;

import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.genotype.mappers.NumDistinctAllelesLabelMapper;
import org.campagnelab.dl.genotype.mappers.SingleBaseFeatureMapperV1;
import org.campagnelab.dl.genotype.mappers.SingleBaseLabelMapperV1;
import org.campagnelab.dl.genotype.tools.SegmentTrainingArguments;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.lossfunctions.ILossFunction;

import java.io.IOException;
import java.util.Map;
import java.util.Properties;
import java.util.function.Function;

/** A domain descriptor for the task of predicting genotypes for segments/blocks of the genome.
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

        if (cachedFeatureMappers.containsKey(inputName)) {
            return cachedFeatureMappers.get(inputName);
        } else {
            switch (inputName) {
                case "input":
                    return new SingleBaseFeatureMapperV1(0);
                default:
                    throw new RuntimeException("Unsupported input name: " + inputName);
            }
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
                    return new SingleBaseLabelMapperV1(0);
                default:
                    throw new RuntimeException("Unsupported output name: " + outputName);
            }
        }
    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {
        return null;
    }

    @Override
    public Function<String, ? extends Iterable<SegmentInformationRecords.SegmentInformation>> getRecordIterable() {
        return inputFilename -> {
            try {
                return new SegmentBaseInformationReader(inputFilename);
            } catch (IOException e) {
                throw new RuntimeException("Unable to read records from " + inputFilename, e);
            }
        };
    }

    @Override
    public ComputationGraphAssembler getComputationalGraph() {
        return null;
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
        // same as input size.
        return getNumInputs(inputName);
    }

    @Override
    public int[] getNumMaskOutputs(String outputName) {
        // same as input size.
        return getNumOutputs(outputName);
    }
    @Override
    public int getNumHiddenNodes(String componentName) {
        return 1024;
    }

    @Override
    public ILossFunction getOutputLoss(String outputName) {
        return null;
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        return 0;
    }
}
