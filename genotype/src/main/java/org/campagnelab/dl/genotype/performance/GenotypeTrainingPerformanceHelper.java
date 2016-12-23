package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.PredictWithModel;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

/**

 */
public class GenotypeTrainingPerformanceHelper extends PredictWithModel<BaseInformationRecords.BaseInformation> {

    protected StatsAccumulator accumulator;


    public GenotypeTrainingPerformanceHelper(DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor) {
        super(domainDescriptor);
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graph,
                                    Predicate<Integer> stopIfTrue) {
        return estimateWithGraph(iterator, graph, stopIfTrue, genotypePrediction -> {
        }, score -> {
        });
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graph,
                                    Predicate<Integer> stopIfTrue,
                                    Consumer<GenotypePrediction > observer, Consumer<Double> scoreObserver) {
        iterator.reset();
        accumulator = new StatsAccumulator();
        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;

        List<Prediction> predictions = new ArrayList<>();
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            double dsScore = graph.score();
            if (dsScore == dsScore) {
                scoreObserver.accept(dsScore);
            }
            INDArray[] trueLabels = next.getLabels();

            int numExamples = next.getFeatures(0).size(0);
            for (int predictionIndex = 0; predictionIndex < numExamples; predictionIndex++) {
                predictions.clear();
                for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {

                    if (interpretors[outputIndex] != null) {
                        Prediction prediction = interpretors[outputIndex].interpret(
                                trueLabels[outputIndex],
                                outputs[outputIndex],
                                predictionIndex);
                        prediction.outputIndex = outputIndex;
                        prediction.index = index;
                        predictions.add(prediction);
                    }
                }
                GenotypePrediction gp = (GenotypePrediction) domainDescriptor.aggregatePredictions(predictions);
                accumulator.observe(gp);
                observer.accept(gp);
                if (stopIfTrue.test(nProcessed)) {
                    break;
                }
                nProcessed += 1;
            }
        }

        return accumulator.createOutputStatistics()[StatsAccumulator.F1_INDEX];
    }

    public double[] getMetricValues(String... metrics) {
        return accumulator.createOutputStatistics(metrics);
    }
}

