package org.campagnelab.dl.somatic.learning.domains;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.performance.AUCHelper;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.somatic.learning.SomaticTrainer;
import org.campagnelab.dl.somatic.learning.SomaticTrainingArguments;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.somatic.mappers.IsBaseMutatedMapper;
import org.campagnelab.dl.somatic.mappers.IsSomaticMutationMapper;
import org.campagnelab.dl.somatic.mappers.SomaticFrequencyLabelMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.ILossFunction;
import org.nd4j.linalg.lossfunctions.impl.LossMCXENT;
import org.nd4j.linalg.lossfunctions.impl.LossMSE;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SomaticMutationDomainDescriptor extends DomainDescriptor<BaseInformationRecords.BaseInformation> {


    private float reductionRate;
    private int genomicContextSize;
    private int indelSequenceLength;
    private float modelCapacity;
    private int ploidy;

    public SomaticMutationDomainDescriptor(SomaticTrainingArguments arguments) {
        this.arguments = arguments;
        initializeArchitecture(arguments.architectureClassname);
    }

    /**
     * Use this method to create a domain descriptor for a trained model.
     *
     * @param modelPath Path where the model is stored.
     */
    public SomaticMutationDomainDescriptor(String modelPath) {

        this.arguments = new SomaticTrainingArguments();
        super.loadProperties(modelPath);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
        this.ploidy=Integer.parseInt(modelProperties.getProperty("genotypes.ploidy"));
        this.genomicContextSize=Integer.parseInt(modelProperties.getProperty("stats.genomicContextSize.min"));
        this.indelSequenceLength=Integer.parseInt(modelProperties.getProperty("indelSequenceLength"));
        this.modelCapacity=Float.parseFloat(modelProperties.getProperty("modelCapacity"));
        this.reductionRate=Float.parseFloat(modelProperties.getProperty("reductionRate"));
        initializeArchitecture();

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

            genomicContextSize = args().genomicContextLength;
            modelProperties.setProperty("stats.genomicContextSize.min", Integer.toString(args().genomicContextLength));
            modelProperties.setProperty("stats.genomicContextSize.max", Integer.toString(args().genomicContextLength));

            indelSequenceLength = args().indelSequenceLength;
            modelProperties.setProperty("indelSequenceLength", Integer.toString(args().indelSequenceLength));

            modelCapacity = args().modelCapacity;
            modelProperties.setProperty("modelCapacity", Float.toString(args().modelCapacity));

            reductionRate = args().reductionRate;
            modelProperties.setProperty("reductionRate", Float.toString(args().reductionRate));

            modelProperties.setProperty("labelSmoothing.epsilon", Double.toString(args().labelSmoothingEpsilon));

            ploidy = args().ploidy;
            modelProperties.setProperty("genotypes.ploidy", Integer.toString(args().ploidy));
            modelProperties.setProperty("genotypes.ploidy", Integer.toString(args().ploidy));
        }
    }

    /**
     * Use this method to create a domain before training. The supplied properties provide
     * featureMapper and labelMapper information (domainProperties), and the sbiProperties
     * provide statistic observed on the training set (or part of it).
     *
     * @param domainProperties Properties describing the domain. Must describe feature and label mappers.
     * @param sbiProperties    Properties describing statistics, used to configure feature mappers.
     */
    public SomaticMutationDomainDescriptor(Properties domainProperties, Properties sbiProperties) {

        this.arguments = new SomaticTrainingArguments();
        super.loadProperties(domainProperties, sbiProperties);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
        initializeArchitecture();
    }

    private SomaticTrainingArguments arguments;

    private SomaticTrainingArguments args() {
        return arguments;
    }

    Map<String, FeatureMapper> cachedFeatureMappers = new HashMap<>();
    public  FeatureMapper configureFeatureMapper(String featureMapperClassname, boolean isTrio, String[] trainingSets) throws IOException {


        try {
            Class clazz = Class.forName(featureMapperClassname + (isTrio ? "Trio" : ""));
            final FeatureMapper featureMapper = (FeatureMapper) clazz.newInstance();
            if (featureMapper instanceof ConfigurableFeatureMapper) {
                ConfigurableFeatureMapper cmapper = (ConfigurableFeatureMapper) featureMapper;
                final Properties sbiProperties = SomaticTrainer.getReaderProperties(trainingSets[0]);
                decorateProperties(sbiProperties);
                cmapper.configure(sbiProperties);
            }
            return featureMapper;
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        if (cachedFeatureMappers.containsKey(inputName)) {
            return cachedFeatureMappers.get(inputName);
        }
        FeatureMapper result;
        if (args().featureMapperClassname != null) {
            assert "input".equals(inputName) : "Only one input supported by this domain.";
            try {

                result = configureFeatureMapper(args().featureMapperClassname, (args()).isTrio,
                        args().getTrainingSets());

            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } else {
            try {
                FeatureMapper fMapper = (FeatureMapper) Class.forName(domainProperties.getProperty("input.featureMapper")).newInstance();
                if (fMapper instanceof ConfigurableFeatureMapper) {
                    ConfigurableFeatureMapper cfmapper = (ConfigurableFeatureMapper) fMapper;
                    cfmapper.configure(modelProperties);
                }
                result = fMapper;
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        cachedFeatureMappers.put(inputName, result);
        return result;

    }
    @Override
    public String produceCacheUniqueId(int miniBatchSize) {
        String id = super.produceCacheUniqueId(miniBatchSize);
        int domainHashcode = id.hashCode();
        domainHashcode ^= genomicContextSize;
        domainHashcode ^= Float.hashCode(args().labelSmoothingEpsilon);
        domainHashcode ^= Integer.hashCode(args().indelSequenceLength);
        domainHashcode ^= Integer.hashCode(args().ploidy);
        return Integer.toHexString(domainHashcode);
    }

    @Override
    public LabelMapper getLabelMapper(String outputName) {

        switch (outputName) {
            case "isMutated":
                return new IsSomaticMutationMapper();
            case "isBaseMutated":
                return new IsBaseMutatedMapper(ploidy);
            case "somaticFrequency":
                return new SomaticFrequencyLabelMapper();
            default:
                throw new IllegalArgumentException("output name is not recognized: " + outputName);
        }
    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {
        switch (outputName) {
            case "isMutated":
                return new IsSomaticMutationInterpreter();
            case "isBaseMutated":
                return new IsBaseMutatedInterpreter();
            case "somaticFrequency":
                return new SomaticFrequencyInterpreter();
            default:
                throw new IllegalArgumentException("output name is not recognized: " + outputName);
        }
    }

    @Override
    public Function<String, ? extends Iterable<BaseInformationRecords.BaseInformation>> getRecordIterable() {
        return inputFilename -> {
            try {
                return new SequenceBaseInformationReader(inputFilename);
            } catch (IOException e) {
                throw new RuntimeException("Unable to read records from " + inputFilename, e);
            }
        };
    }

    @Override
    public PerformanceMetricDescriptor<BaseInformationRecords.BaseInformation> performanceDescritor() {
        return new PerformanceMetricDescriptor<BaseInformationRecords.BaseInformation>(this) {
            @Override
            public String[] performanceMetrics() {
                return new String[]{"AUC", "score"};
            }

            @Override
            public boolean largerValueIsBetterPerformance(String metricName) {
                switch (metricName) {
                    case "AUC":
                        return true;
                    case "score":
                        return false;
                    default:
                        throw new IllegalArgumentException("metric not recognized.");
                }
            }

            @Override
            public double estimateMetric(ComputationGraph graph, String metricName, MultiDataSetIterator dataSetIterator, long scoreN) {
                switch (metricName) {
                    case "AUC":
                        AUCHelper helper = new AUCHelper();
                        return helper.estimateWithGraph(dataSetIterator, graph, args().numValidation, prediction -> {
                                },
                                index -> index > scoreN,
                            /* first output represents probability of mutation */ 0,
                                hasOutput("isBaseMutated") ? new IsBaseMutatedInterpreter() :
                                        new IsSomaticMutationInterpreter());
                    default:
                        return estimateScore(graph, metricName, dataSetIterator, scoreN);
                }
            }

            @Override
            public String earlyStoppingMetric() {
                return args().earlyStoppingMeasureName;
            }
        };


    }


    /**
     * Initialize the net architecture using the name in the domain properties.
     */
    @Override
    protected void initializeArchitecture() {
        String netArchitectureClassname = domainProperties.getProperty("net.architecture.classname");
        if (netArchitectureClassname == null) {
            netArchitectureClassname = SixDenseLayersNarrower2.class.getCanonicalName();
        }
        initializeArchitecture(netArchitectureClassname);
    }

    @Override
    public ComputationGraphAssembler getComputationalGraph() {
        computationGraphAssembler.setArguments(arguments);
        return computationGraphAssembler;
    }

    @Override
    public int[] getNumInputs(String inputName) {
        return new int[]{getFeatureMapper(inputName).numberOfFeatures()};
    }

    @Override
    public int[] getNumOutputs(String outputName) {
        return new int[]{getLabelMapper(outputName).numberOfLabels()};
    }

    // TODO: SomaticMutationDomainDescriptor shouldn't need these methods, but make sure
    @Override
    public int[] getNumMaskInputs(String inputName) {
        return new int[]{getFeatureMapper(inputName).numberOfFeatures()};
    }

    @Override
    public int[] getNumMaskOutputs(String outputName) {
        return new int[]{getLabelMapper(outputName).numberOfLabels()};
    }

    @Override
    public int getNumHiddenNodes(String componentName) {
        return (int) (getNumInputs("input")[0] *modelCapacity);
    }

    @Override
    public ILossFunction getOutputLoss(String outputName) {
        switch (outputName) {
            case "isMutated":
            case "isBaseMutated":
                return new LossMCXENT();
            case "somaticFrequency":
                return new LossMSE();
            default:
                throw new IllegalArgumentException("Output name is not recognized");
        }
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        BaseInformationConcatIterator it = null;
        try {
            List<BaseInformationIterator> list = Arrays.asList(recordFiles).stream().map(filename -> {
                try {

                    return new BaseInformationIterator(filename, 128, featureMappers()[0], getLabelMapper("isMutated"));
                } catch (IOException e) {
                    throw new RuntimeException("Unable to estimate number of records for filename " + filename);
                }
            }).collect(Collectors.toList());
            it = new BaseInformationConcatIterator(list, 128, featureMappers()[0], getLabelMapper("isMutated"));
            return it.totalExamples();

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            IOUtils.closeQuietly(it);
        }
        return 0;
    }
    @Override
    public void putProperties(Properties props) {
        super.putProperties(props);
        decorateProperties(props);
    }
}