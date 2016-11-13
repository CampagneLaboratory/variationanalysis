package org.campagnelab.dl.varanalysis.learning.domains;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.model.utils.ConfigurableFeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.IsSomaticMutationMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SomaticFrequencyLabelMapper;
import org.campagnelab.dl.varanalysis.learning.SomaticTrainingArguments;
import org.campagnelab.dl.varanalysis.learning.TrainSomaticModel;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.campagnelab.dl.varanalysis.learning.architecture.graphs.SixDenseLayersNarrower2;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationConcatIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AUCHelper;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SomaticMutationDomainDescriptor extends DomainDescriptor<BaseInformationRecords.BaseInformation> {


    public SomaticMutationDomainDescriptor(SomaticTrainingArguments arguments) {
        this.arguments = arguments;
    }

    public SomaticMutationDomainDescriptor(String modelPath) {

        this.arguments = new SomaticTrainingArguments();
        super.loadProperties(modelPath);
        // force loading the feature mappers from properties.
        args().featureMapperClassname = null;
    }

    private SomaticTrainingArguments arguments;

    private SomaticTrainingArguments args() {
        return arguments;
    }

    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        if (args().featureMapperClassname != null) {
            assert "input".equals(inputName) : "Only one input supported by this domain.";

            try {
                return TrainSomaticModel.configureFeatureMapper(args().featureMapperClassname, (args()).isTrio,
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


    @Override
    public LabelMapper getLabelMapper(String outputName) {

        switch (outputName) {
            case "isMutated":
                return new IsSomaticMutationMapper();
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
        return new SixDenseLayersNarrower2();
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
        return getNumInputs("input")[0] * 4;
    }

    @Override
    public LossFunctions.LossFunction getOutputLoss(String outputName) {
        switch (outputName) {
            case "isMutated":
                return LossFunctions.LossFunction.MCXENT;
            case "somaticFrequency":
                return LossFunctions.LossFunction.MSE;
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

}