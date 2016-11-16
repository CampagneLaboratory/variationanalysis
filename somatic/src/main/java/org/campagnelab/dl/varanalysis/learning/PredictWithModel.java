package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.varanalysis.models.ModelOutputHelper;
import org.campagnelab.dl.varanalysis.learning.domains.DomainDescriptor;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.Prediction;
import org.campagnelab.dl.varanalysis.learning.domains.predictions.PredictionInterpreter;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

/**
 * Helper class to predict with a model and obtain interpreted predictions.
 */
public class PredictWithModel<RecordType> {
    DomainDescriptor<RecordType> domainDescriptor;
    ModelOutputHelper outputHelper;
    private PredictionInterpreter[] interpretors;

    public PredictWithModel(DomainDescriptor<RecordType> domainDescriptor) {
        this.domainDescriptor = domainDescriptor;
        outputHelper = new ModelOutputHelper();
        int outputIndex = 0;

        String[] outputNames = domainDescriptor.getComputationalGraph().getOutputNames();
        interpretors = new PredictionInterpreter[outputNames.length];
        for (String outputName : outputNames) {
            interpretors[outputIndex++] = domainDescriptor.getPredictionInterpreter(outputName);
        }
    }

    public void makePredictions(Iterator<RecordType> iterator,
                                Model model,
                                Consumer<List<Prediction>> doForEachPrediction,
                                Predicate<Integer> stopIfTrue) {
        makePredictions(iterator, model,
                recordType -> {
        },
                doForEachPrediction,
                stopIfTrue);
    }

    public void makePredictions(Iterator<RecordType> iterator,
                                Model model,
                                Consumer<RecordType> observeRecord,
                                Consumer<List<Prediction>> doForEachPrediction,
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
            doForEachPrediction.accept(predictions);
            index++;
            if (stopIfTrue.test(index)) {
                break;
            }
        }
    }


}
