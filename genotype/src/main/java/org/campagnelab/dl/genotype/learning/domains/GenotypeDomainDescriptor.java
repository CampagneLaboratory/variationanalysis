package org.campagnelab.dl.genotype.learning.domains;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.framework.architecture.graphs.ComputationalGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.performance.AUCHelper;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.genotype.learning.GenotypeTrainingArguments;
import org.campagnelab.dl.genotype.learning.TrainGenotypeModelS;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2;
import org.campagnelab.dl.genotype.learning.domains.predictions.GenotypeInterpreter;
import org.campagnelab.dl.genotype.mappers.GenotypeLabelsMapper;
import org.campagnelab.dl.somatic.learning.SomaticTrainingArguments;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2;
import org.campagnelab.dl.somatic.learning.domains.SomaticFrequencyInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.somatic.mappers.IsSomaticMutationMapper;
import org.campagnelab.dl.somatic.mappers.SomaticFrequencyLabelMapper;
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

public class GenotypeDomainDescriptor extends DomainDescriptor<BaseInformationRecords.BaseInformation> {


    public GenotypeDomainDescriptor(GenotypeTrainingArguments arguments) {
        this.arguments = arguments;
    }

    /**
     * Use this method to create a domain descriptor for a trained model.
     * @param modelPath Path where the model is stored.
     */
    public GenotypeDomainDescriptor(String modelPath) {

        this.arguments = new GenotypeTrainingArguments();
        super.loadProperties(modelPath);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
    }

    /**
     * Use this methdo to create a domain before training. The supplied properties provide
     * featureMapper and labelMapper information (domainProperties), and the sbiProperties
     * provide statistic observed on the training set (or part of it).
     * @param domainProperties Properties describing the domain. Must describe feature and label mappers.
     * @param sbiProperties Properties describing statistics, used to configure feature mappers.
     */
    public GenotypeDomainDescriptor(Properties domainProperties, Properties sbiProperties) {

        this.arguments = new GenotypeTrainingArguments();
        super.loadProperties(domainProperties,sbiProperties);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
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


        switch (outputName) {
            case "A":
                return new GenotypeLabelsMapper(0);
            case "T":
                return new GenotypeLabelsMapper(1);
            case "C":
                return new GenotypeLabelsMapper(2);
            case "G":
                return new GenotypeLabelsMapper(3);
            case "N":
                return new GenotypeLabelsMapper(4);
            case "I1":
                return new GenotypeLabelsMapper(5);
            case "I2":
                return new GenotypeLabelsMapper(6);
            case "I3":
                return new GenotypeLabelsMapper(7);
            case "I4":
                return new GenotypeLabelsMapper(8);
            case "I5":
                return new GenotypeLabelsMapper(9);
            case "genotype":
                //handle this case for the properties file description
                return new GenotypeLabelsMapper(0);

            default:
                throw new IllegalArgumentException("output name is not recognized: " + outputName);
        }

    }

    @Override
    public PredictionInterpreter getPredictionInterpreter(String outputName) {

        switch (outputName) {
            case "A":
                return new GenotypeInterpreter(0);
            case "T":
                return new GenotypeInterpreter(1);
            case "C":
                return new GenotypeInterpreter(2);
            case "G":
                return new GenotypeInterpreter(3);
            case "N":
                return new GenotypeInterpreter(4);
            case "I1":
                return new GenotypeInterpreter(5);
            case "I2":
                return new GenotypeInterpreter(6);
            case "I3":
                return new GenotypeInterpreter(7);
            case "I4":
                return new GenotypeInterpreter(8);
            case "I5":
                return new GenotypeInterpreter(9);
            case "genotype":
                //handle this case for the properties file description
                return new GenotypeInterpreter(0);

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
                            /* first output represents probability of mutation */ 0);
                    default:
                        return estimateScore(graph, metricName, dataSetIterator, scoreN);
                }
            }

            @Override
            public String earlyStoppingMetric() {
                return "AUC";
            }
        };


    }

    @Override
    public ComputationalGraphAssembler getComputationalGraph() {
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
        return LossFunctions.LossFunction.XENT;
    }

    @Override
    public long getNumRecords(String[] recordFiles) {
        BaseInformationConcatIterator it = null;
        try {
            List<BaseInformationIterator> list = Arrays.asList(recordFiles).stream().map(filename -> {
                try {

                    return new BaseInformationIterator(filename, 128, featureMappers()[0], getLabelMapper("genotype"));
                } catch (IOException e) {
                    throw new RuntimeException("Unable to estimate number of records for filename " + filename);
                }
            }).collect(Collectors.toList());
            it = new BaseInformationConcatIterator(list, 128, featureMappers()[0], getLabelMapper("genotype"));
            return it.totalExamples();

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            IOUtils.closeQuietly(it);
        }
        return 0;
    }

}