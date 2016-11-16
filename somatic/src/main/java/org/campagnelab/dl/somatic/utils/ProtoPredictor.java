package org.campagnelab.dl.somatic.utils;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
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

    public static final int POSITIVE_STRAND = 0;
    public static final int NEGATIVE_STRAND = 1;
    public static final int POSITIVE_PROBABILITY_INDEX = 0;
    public static final int NEGATIVE_PROBABILITY_INDEX = 1;
    private Model model;
    private FeatureMapper mapper;
    private ModelOutputHelper outputHelper;

    public ProtoPredictor(Model model, FeatureMapper mapper) {
        this.model = model;
        this.mapper = mapper;
        this.outputHelper = new ModelOutputHelper();
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



    public Prediction mutPrediction(BaseInformationRecords.BaseInformation record) {
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
            arrayPredicted = outputHelper.getOutput(0);
            float[] probabilities = arrayPredicted.getRow(0).data().asFloat();
            prediction.set(probabilities[POSITIVE_PROBABILITY_INDEX],
                    probabilities[NEGATIVE_PROBABILITY_INDEX]);
            prediction.setPredictedSomaticFrequency( outputHelper.getOutput(1).getFloat(0));
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


