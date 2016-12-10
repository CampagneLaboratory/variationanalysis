package org.campagnelab.dl.somatic.utils;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.campagnelab.dl.somatic.learning.domains.SomaticFrequencyInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsMutatedPrediction;
import org.campagnelab.dl.somatic.learning.domains.predictions.IsSomaticMutationInterpreter;
import org.campagnelab.dl.somatic.learning.domains.predictions.SomaticFrequencyPrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.List;

/**
 * Created by rct66 on 6/23/16.
 */
public class ProtoPredictor {

    public static final int POSITIVE_PROBABILITY_INDEX = 0;
    public static final int NEGATIVE_PROBABILITY_INDEX = 1;
    private static final int PROBABILITY_OUTPUT_INDEX = 0;
    private static final int SOMATIC_FREQUENCY_INDEX = 1;
    private Model model;
    private FeatureMapper mapper;
    private ModelOutputHelper outputHelper;
    private DomainDescriptor domainDescriptor;

    public ProtoPredictor(DomainDescriptor domainDescriptor, Model model, FeatureMapper mapper) {
        this.model = model;
        this.mapper = mapper;
        this.outputHelper = new ModelOutputHelper();
        this.domainDescriptor = domainDescriptor;
        if (domainDescriptor != null) {
            if (domainDescriptor.hasOutput("isMutated")) {
                isSomatic = (IsSomaticMutationInterpreter) domainDescriptor.getPredictionInterpreter("isMutated");
            }
            if (domainDescriptor.hasOutput("isBaseMutated")) {
                isSomatic = (IsSomaticMutationInterpreter) domainDescriptor.getPredictionInterpreter("isBaseMutated");
            }
            if (domainDescriptor.hasOutput("somaticFrequency")) {
                somaticFrequency = (SomaticFrequencyInterpreter) domainDescriptor.getPredictionInterpreter("somaticFrequency");
            }
        } else {
            // for backward compatibility with models not trained with the framework:
            isSomatic = new IsSomaticMutationInterpreter();
            somaticFrequency = new SomaticFrequencyInterpreter();
        }

    }

    public static List<Integer> expandFreq(List<BaseInformationRecords.NumberWithFrequency> freqList) {
        int capacity = 0;
        for (BaseInformationRecords.NumberWithFrequency freq : freqList) {
            capacity += freq.getFrequency();
        }
        IntArrayList expanded = new IntArrayList(capacity);
        for (BaseInformationRecords.NumberWithFrequency freq : freqList) {
            for (int i = 0; i < freq.getFrequency(); i++) {
                expanded.add(freq.getNumber());
            }
        }

        return expanded;
    }

    PredictionInterpreter<BaseInformationRecords.BaseInformation, IsMutatedPrediction> isSomatic = null;
    SomaticFrequencyInterpreter somaticFrequency = null;

    public Prediction mutPrediction(BaseInformationRecords.BaseInformation record) {
        assert model != null : "Model cannot be null";
        assert isSomatic != null : "isSomatic interpreter must not be null";
        assert somaticFrequency != null : "somaticFrequency interpreter must not be null";
        INDArray arrayPredicted = null;
        Prediction prediction = new Prediction();

        if (model instanceof MultiLayerNetwork) {
            INDArray testFeatures = Nd4j.zeros(1, mapper.numberOfFeatures());
            mapper.prepareToNormalize(record, 0);
            mapper.mapFeatures(record, testFeatures, 0);
            arrayPredicted = ((MultiLayerNetwork) model).output(testFeatures, false);
            float[] probabilities = arrayPredicted.getRow(0).data().asFloat();
            prediction.set(probabilities[POSITIVE_PROBABILITY_INDEX],
                    probabilities[NEGATIVE_PROBABILITY_INDEX]);

        } else if (model instanceof ComputationGraph) {
            outputHelper.predictForNextRecord(model, record, mapper);

            IsMutatedPrediction isSomaticPrediction = isSomatic.interpret(record, outputHelper.getOutput(PROBABILITY_OUTPUT_INDEX));
            prediction.set((float) isSomaticPrediction.predictedLabelYes,
                    (float) isSomaticPrediction.predictedLabelNo);
            SomaticFrequencyPrediction somaticFrequencyPrediction = somaticFrequency.interpret(record,
                    outputHelper.getOutput(SOMATIC_FREQUENCY_INDEX));
            prediction.setPredictedSomaticFrequency(somaticFrequencyPrediction.predictedValue);
        }

        return prediction;
    }

    public Prediction getNullPrediction() {
        return new Prediction(0, 0);
    }

    public class Prediction {
        public boolean clas;
        public float posProb;
        public float negProb;
        public float predictedSomaticFrequency;
        public boolean hasSomaticFrequency;

        public boolean isMutated() {
            return posProb > negProb;
        }

        public Prediction(float posProb, float negProb) {
            set(posProb, negProb);
        }

        public Prediction() {
        }

        public void set(float posProb, float negProb) {
            this.posProb = posProb;
            this.negProb = negProb;
            this.clas = (posProb > 0.5);
        }

        public void setPredictedSomaticFrequency(float value) {
            hasSomaticFrequency = true;
            this.predictedSomaticFrequency = value;
        }

        public boolean isCorrect(boolean isPositiveTrueLabel) {
            return isMutated() && isPositiveTrueLabel || (!isMutated() && !isPositiveTrueLabel);
        }
    }
}


