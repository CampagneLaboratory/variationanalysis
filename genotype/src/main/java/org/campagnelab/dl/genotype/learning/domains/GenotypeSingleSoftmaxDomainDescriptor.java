package org.campagnelab.dl.genotype.learning.domains;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2;
import org.campagnelab.dl.genotype.learning.domains.predictions.HomozygousInterpreter;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypeInterpreter;
import org.campagnelab.dl.genotype.mappers.GenotypeFeatureMapper;
import org.campagnelab.dl.genotype.mappers.GenotypeLabelsMapper;
import org.campagnelab.dl.genotype.mappers.HomozygousLabelsMapper;
import org.campagnelab.dl.genotype.performance.AccuracyHelper;
import org.campagnelab.dl.genotype.performance.AlleleAccuracyHelper;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.function.Function;
import java.util.stream.Collectors;

public class GenotypeSingleSoftmaxDomainDescriptor extends DomainDescriptor<BaseInformationRecords.BaseInformation> {


    public GenotypeSingleSoftmaxDomainDescriptor(GenotypeTrainingArguments arguments) {
        this.arguments = arguments;
        initializeArchitecture(arguments.architectureClassname);
    }

    /**
     * Use this method to create a domain descriptor for a trained model.
     *
     * @param modelPath Path where the model is stored.
     */
    public GenotypeSingleSoftmaxDomainDescriptor(String modelPath) {

        this.arguments = new GenotypeTrainingArguments();
        super.loadProperties(modelPath);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
        initializeArchitecture();
    }

    /**
     * Use this methdo to create a domain before training. The supplied properties provide
     * featureMapper and labelMapper information (domainProperties), and the sbiProperties
     * provide statistic observed on the training set (or part of it).
     *
     * @param domainProperties Properties describing the domain. Must describe feature and label mappers.
     * @param sbiProperties    Properties describing statistics, used to configure feature mappers.
     */
    public GenotypeSingleSoftmaxDomainDescriptor(Properties domainProperties, Properties sbiProperties) {

        this.arguments = new GenotypeTrainingArguments();
        super.loadProperties(domainProperties, sbiProperties);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
        initializeArchitecture();
    }

    private GenotypeTrainingArguments arguments;

    private GenotypeTrainingArguments args() {
        return arguments;
    }

    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        if (args().featureMapperClassname != null) {
            assert "input".equals(inputName) : "Only one input supported by this domain.";

            try {
                return TrainSomaticModel.configureFeatureMapper(args().featureMapperClassname, false,
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
                return fMapper;
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }

    }

    // TODO: SomaticMutationDomainDescriptor shouldn't need these methods, but make sure
    public int[] getNumMaskInputs(String inputName) {
        return new int[]{getFeatureMapper(inputName).numberOfFeatures()};
    }


    @Override
    public int[] getNumMaskOutputs(String outputName) {
        return new int[]{getLabelMapper(outputName).numberOfLabels()};
    }

    @Override
    public LabelMapper getLabelMapper(String outputName) {
        boolean sortCounts=needSortCounts();

        switch (outputName) {
            case "A":
                return new GenotypeLabelsMapper(0,sortCounts);
            case "T":
                return new GenotypeLabelsMapper(1,sortCounts);
            case "C":
                return new GenotypeLabelsMapper(2,sortCounts);
            case "G":
                return new GenotypeLabelsMapper(3,sortCounts);
            case "N":
                return new GenotypeLabelsMapper(4,sortCounts);
            case "I1":
                return new GenotypeLabelsMapper(5,sortCounts);
            case "I2":
                return new GenotypeLabelsMapper(6,sortCounts);
            case "I3":
                return new GenotypeLabelsMapper(7,sortCounts);
            case "I4":
                return new GenotypeLabelsMapper(8,sortCounts);
            case "I5":
                return new GenotypeLabelsMapper(9,sortCounts);
            case "homozygous":
                return new HomozygousLabelsMapper(sortCounts);
            default:
                throw new IllegalArgumentException("output name is not recognized: " + outputName);
        }

    }

    private boolean needSortCounts() {
        return ((GenotypeFeatureMapper)getFeatureMapper("input")).sortCounts;
    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {
        boolean sortCounts=needSortCounts();

        switch (outputName) {
            case "A":
                return new SingleGenotypeInterpreter(0, sortCounts);
            case "T":
                return new SingleGenotypeInterpreter(1, sortCounts);
            case "C":
                return new SingleGenotypeInterpreter(2, sortCounts);
            case "G":
                return new SingleGenotypeInterpreter(3, sortCounts);
            case "N":
                return new SingleGenotypeInterpreter(4, sortCounts);
            case "I1":
                return new SingleGenotypeInterpreter(5, sortCounts);
            case "I2":
                return new SingleGenotypeInterpreter(6, sortCounts);
            case "I3":
                return new SingleGenotypeInterpreter(7, sortCounts);
            case "I4":
                return new SingleGenotypeInterpreter(8, sortCounts);
            case "I5":
                return new SingleGenotypeInterpreter(9, sortCounts);
            //only need one interpreter for each record, it will collect entire genotype into a prediction
            case "homozygous":
                return new HomozygousInterpreter(sortCounts);
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
        return new PerformanceMetricDescriptor<BaseInformationRecords.BaseInformation>() {
            @Override
            public String[] performanceMetrics() {
                return new String[]{"accuracy", "alleleAccuracy", "score"};
            }

            @Override
            public boolean largerValueIsBetterPerformance(String metricName) {
                switch (metricName) {
                    case "accuracy":
                        return true;
                    case "alleleAccuracy":
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
                    case "accuracy":
                        AccuracyHelper helper = new AccuracyHelper();
                        return helper.estimateWithGraph(dataSetIterator, graph,
                                index -> index > scoreN
                            /* first output represents probabilityIsCalled of mutation */);
                    case "alleleAccuracy":
                        AlleleAccuracyHelper alleleHelper = new AlleleAccuracyHelper();
                        return alleleHelper.estimateWithGraph(dataSetIterator, graph,
                                index -> index > scoreN
                            /* first output represents probabilityIsCalled of mutation */);
                    default:
                        return estimateScore(graph, metricName, dataSetIterator, scoreN);
                }
            }

            @Override
            public String earlyStoppingMetric() {
                return "score";
            }
        };


    }

    @Override
    public ComputationGraphAssembler getComputationalGraph() {
        return new GenotypeSixDenseLayersNarrower2();
    }

    @Override
    public int[] getNumInputs(String inputName) {
        return new int[]{getFeatureMapper(inputName).numberOfFeatures()};
    }

    @Override
    public int[] getNumOutputs(String outputName) {
        return new int[]{getLabelMapper(outputName).numberOfLabels()};
    }

    @Override
    public int getNumHiddenNodes(String componentName) {
        return getNumInputs("input")[0] * 1;
    }

    @Override
    public LossFunctions.LossFunction getOutputLoss(String outputName) {
        switch (outputName) {
            case "homozygous":
                return LossFunctions.LossFunction.MCXENT;
            default:
                // any other is an individual genotype output:
                return LossFunctions.LossFunction.XENT;
             }
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        BaseInformationConcatIterator it = null;
        try {
            List<BaseInformationIterator> list = Arrays.asList(recordFiles).stream().map(filename -> {
                try {

                    return new BaseInformationIterator(filename, 128, featureMappers()[0], getLabelMapper("A"));
                } catch (IOException e) {
                    throw new RuntimeException("Unable to estimate number of records for filename " + filename);
                }
            }).collect(Collectors.toList());
            it = new BaseInformationConcatIterator(list, 128, featureMappers()[0], getLabelMapper("A"));
            return it.totalExamples();

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            IOUtils.closeQuietly(it);
        }
        return 0;
    }

}