package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.domains.prediction.RecordPredictions;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

/**
 * Helper class to predict with a model and obtain interpreted predictions.
 */
public class PredictWithModel<RecordType> {
    protected DomainDescriptor<RecordType> domainDescriptor;
    ModelOutputHelper outputHelper;
    protected PredictionInterpreter[] interpretors;

    public PredictWithModel(DomainDescriptor<RecordType> domainDescriptor) {
        this.domainDescriptor = domainDescriptor;
        outputHelper = new ModelOutputHelper<RecordType>(domainDescriptor);
        int outputIndex = 0;

        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();
        interpretors = new PredictionInterpreter[outputNames.length];
        for (String outputName : outputNames) {
            interpretors[outputIndex++] = domainDescriptor.getPredictionInterpreter(outputName);
        }
    }

    public void makePredictions(Iterator<RecordType> iterator,
                                Model model,
                                Consumer<RecordPredictions<RecordType>> doForEachPrediction,
                                Predicate<Integer> stopIfTrue) {
        makePredictions(iterator, model,
                recordType -> {
        },
                doForEachPrediction,
                stopIfTrue);
    }

    public int makePredictions(MultiDataSet dataSet, List<RecordType> records,
                                Model model,
                                Consumer<RecordPredictions<RecordType>> doForEachPrediction,
                                Predicate<Integer> stopIfTrue, int index) {
        assert model instanceof ComputationGraph : "MultiDataSet only work with ComputationGraph";
        ComputationGraph graph=(ComputationGraph)model;
        INDArray[] outputPredictions = graph.output(dataSet.getFeatures());
        List<Prediction> predictions = new ArrayList<>();

        RecordType currentRecord;
        for (int exampleIndex = 0; exampleIndex<records.size(); exampleIndex++) {
            predictions.clear();
            currentRecord=records.get(exampleIndex);
            for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {


                if (interpretors[outputIndex] != null) {
                    Prediction prediction = interpretors[outputIndex].interpret(currentRecord, outputPredictions[outputIndex]);
                    prediction.outputIndex = outputIndex;
                    prediction.index = index;
                    predictions.add(prediction);
                }
            }
            doForEachPrediction.accept(new RecordPredictions<>(currentRecord, predictions));
            index++;
            if (stopIfTrue.test(index)) {
                break;
            }
        }
        return index;

    }

    public void makePredictions(Iterator<RecordType> iterator,
                                Model model,
                                Consumer<RecordType> observeRecord,
                                Consumer<RecordPredictions<RecordType>> doForEachPrediction,
                                Predicate<Integer> stopIfTrue) {
        int index = 0;
        List<Prediction> predictions = new ArrayList<>();
        while (iterator.hasNext()) {

            RecordType currentRecord = iterator.next();
            observeRecord.accept(currentRecord);
            outputHelper.predictForNextRecord(model, currentRecord, domainDescriptor.featureMappers());
            predictions.clear();
            for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {
                INDArray outputPredictions = outputHelper.getOutput(outputIndex);

                if (interpretors[outputIndex] != null) {
                    Prediction prediction = interpretors[outputIndex].interpret(currentRecord, outputPredictions);
                    prediction.outputIndex = outputIndex;
                    prediction.index = index;
                    predictions.add(prediction);
                }
            }

            doForEachPrediction.accept(new RecordPredictions<>(currentRecord, predictions));
            index++;
            if (stopIfTrue.test(index)) {
                break;
            }
        }
    }


}
