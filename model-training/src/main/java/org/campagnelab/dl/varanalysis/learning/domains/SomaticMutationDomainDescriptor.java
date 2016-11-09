package org.campagnelab.dl.varanalysis.learning.domains;

import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.IsSomaticMutationMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.model.utils.mappers.SomaticFrequencyLabelMapper;
import org.campagnelab.dl.varanalysis.learning.DomainDescriptor;
import org.campagnelab.dl.varanalysis.learning.PerformanceMetricDescriptor;
import org.campagnelab.dl.varanalysis.learning.SomaticTrainingArguments;
import org.campagnelab.dl.varanalysis.learning.TrainSomaticModel;
import org.campagnelab.dl.varanalysis.learning.architecture.ComputationalGraphAssembler;
import org.campagnelab.dl.varanalysis.learning.architecture.graphs.SixDenseLayersNarrower2;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AUCHelper;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.lossfunctions.LossFunctions;

import java.io.IOException;
import java.util.function.Function;

public class SomaticMutationDomainDescriptor extends DomainDescriptor<BaseInformationRecords.BaseInformation> {


    public SomaticMutationDomainDescriptor(SomaticTrainingArguments arguments) {
        this.arguments = arguments;
    }

    private SomaticTrainingArguments arguments;

    private SomaticTrainingArguments args() {
        return arguments;
    }

    @Override
    public FeatureMapper getFeatureMapper(String inputName) {
        assert "input".equals(inputName) : "Only one input supported by this domain.";
        FeatureMapper featureMapper = null;
        try {
            featureMapper = TrainSomaticModel.configureFeatureMapper(args().featureMapperClassname, (args()).isTrio,
                    args().getTrainingSets());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return featureMapper;
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
                        return helper.estimate(dataSetIterator, graph, args().numValidation, prediction -> {
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
            SequenceBaseInformationReader reader = null;
            try {
                reader = new SequenceBaseInformationReader(recordFiles[0]);
                return reader.getTotalRecords();
            } catch (IOException e) {
                e.printStackTrace();
            } finally {
                IOUtils.closeQuietly(reader);
            }
            return 0;
        }

}