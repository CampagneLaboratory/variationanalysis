package org.campagnelab.dl.framework.domains.prediction;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelOutputHelper;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

public abstract class PredictWith<RecordType> {

    protected DomainDescriptor<RecordType> domainDescriptor;
    ModelOutputHelper outputHelper;
    protected PredictionInterpreter[] interpretors;
    protected Model model;

    public PredictWith(DomainDescriptor<RecordType> domainDescriptor) {
        this.domainDescriptor = domainDescriptor;
        outputHelper = new ModelOutputHelper<RecordType>(domainDescriptor);
        int outputIndex = 0;

        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();
        interpretors = new PredictionInterpreter[outputNames.length];
        for (String outputName : outputNames) {
            interpretors[outputIndex++] = domainDescriptor.getPredictionInterpreter(outputName);
        }
    }

    public abstract INDArray[] getModelOutputs(MultiDataSet dataSet, int batchSize);

    public int makePredictions(MultiDataSet dataSet, List<RecordType> records,
                               Consumer<RecordPredictions<RecordType>> doForEachPrediction,
                               Predicate<Integer> stopIfTrue, int index) {

        int batchSize=records.size();
        INDArray[] outputPredictions = getModelOutputs(dataSet,batchSize);
        List<Prediction> predictions = new ArrayList<>();

        RecordType currentRecord;
        for (int exampleIndex = 0; exampleIndex<records.size(); exampleIndex++) {
            predictions.clear();
            currentRecord=records.get(exampleIndex);
            for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {


                if (interpretors[outputIndex] != null) {
                    Prediction prediction = interpretors[outputIndex].interpret(currentRecord,
                            outputPredictions[outputIndex].slice(exampleIndex));
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
                                Consumer<RecordType> observeRecord,
                                Consumer<RecordPredictions<RecordType>> doForEachPrediction,
                                Predicate<Integer> stopIfTrue) {
        int index = 0;
        List<Prediction> predictions = new ArrayList<>();
        while (iterator.hasNext()) {

            RecordType currentRecord = iterator.next();
            observeRecord.accept(currentRecord);
            outputHelper.predictForNextRecord(model, currentRecord, domainDescriptor.featureMappers(true));
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
