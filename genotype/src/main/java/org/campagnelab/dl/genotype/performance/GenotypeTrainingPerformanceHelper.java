package org.campagnelab.dl.genotype.performance;

import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.PredictWithModel;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.api.Model;
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


    public GenotypeTrainingPerformanceHelper(DomainDescriptor<BaseInformationRecords.BaseInformation> domainDescriptor, Model graph) {
        super(domainDescriptor,graph);
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graph,
                                    Predicate<Integer> stopIfTrue) {
        return estimateWithGraph(iterator, graph, stopIfTrue, genotypePrediction -> {
        }, score -> {
        });
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graphUnused,
                                    Predicate<Integer> stopIfTrue,
                                    Consumer<GenotypePrediction> observer, Consumer<Double> scoreObserver) {
        iterator.reset();
        accumulator = new StatsAccumulator();
        accumulator.initializeStats();
        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;
        ComputationGraph graph= (ComputationGraph) model;
        List<Prediction> predictions = new ArrayList<>();
        while (iterator.hasNext()) {
            MultiDataSet next = iterator.next();
            INDArray[] outputs = graph.output(next.getFeatures());
            double dsScore = graph.score(next);
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
                GenotypePrediction gp = (GenotypePrediction) domainDescriptor.aggregatePredictions(null/**null in training phase */, predictions);
                // obtain the reference base as an int (e.g., 0 or 1), to match the format of
                // genotypes obtained during training from the cache:
                String referenceBase = Integer.toString(gp.referenceGobyIndex);
                accumulator.observe(gp, gp.isVariant(), GenotypeHelper.isVariant(gp.predictedGenotype, referenceBase));
                observer.accept(gp);
                if (stopIfTrue.test(nProcessed)) {
                    break;
                }
                nProcessed += 1;
                index++;
            }
        }

        return accumulator.createOutputStatistics()[StatsAccumulator.F1_INDEX];
    }

    public double[] getMetricValues(String... metrics) {
        return accumulator.createOutputStatistics(metrics);
    }
}

