package org.campagnelab.dl.genotype.performance;

import it.unimi.dsi.fastutil.doubles.DoubleList;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.PredictWithModel;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

/**

 */
public class SegmentTrainingPerformanceHelper extends PredictWithModel<SegmentInformationRecords.SegmentInformation> {

    protected SegmentStatsAccumulator accumulator;
    double validationScore;
    int n;

    public SegmentTrainingPerformanceHelper(DomainDescriptor<SegmentInformationRecords.SegmentInformation> domainDescriptor) {
        super(domainDescriptor,null);
        accumulator = new SegmentStatsAccumulator(new AccuracySegmentAccumulator(), new IndelAccuracySegmentAccumulator());
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graph,
                                    Predicate<Integer> stopIfTrue) {
        validationScore = 0;
        n=0;
        return estimateWithGraph(iterator, graph, stopIfTrue, genotypePrediction -> {
        }, score -> {
            validationScore += score;
            n++;
        });
    }

    public double estimateWithGraph(MultiDataSetIterator iterator,
                                    ComputationGraph graph,
                                    Predicate<Integer> stopIfTrue,
                                    Consumer<SegmentPrediction> observer, Consumer<Double> scoreObserver) {
        iterator.reset();
        accumulator.initializeStats();
        int index = 0;
        int nProcessed = 0;
        int nCorrect = 0;

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
            for (int exampleIndex = 0; exampleIndex < numExamples; exampleIndex++) {
                predictions.clear();
                for (int outputIndex = 0; outputIndex < domainDescriptor.getNumModelOutputs(); outputIndex++) {

                    if (interpretors[outputIndex] != null) {
                        Prediction prediction = interpretors[outputIndex].interpret(
                                trueLabels[outputIndex].slice(exampleIndex),
                                outputs[outputIndex].slice(exampleIndex),
                                exampleIndex);
                        prediction.outputIndex = outputIndex;
                        prediction.index = index;
                        predictions.add(prediction);
                    }
                }
                SegmentPrediction sp = (SegmentPrediction) domainDescriptor.aggregatePredictions(null/**null in training phase */, predictions);
                accumulator.observe(sp);

                if (stopIfTrue.test(nProcessed)) {
                    break;
                }
                nProcessed += 1;
            }
        }

        return accumulator.estimates().get(0);
    }

    public double[] getMetricValues() {
        DoubleList estimates = accumulator.estimates();
        double[] result = new double[estimates.size() + 1];
        // put the validation score in first position in the result:
        result[0] = validationScore / (double) n;
        for (int i=1;i<estimates.size()+1;i++) {
            result[i]=estimates.get(i-1);
        }
       // Array.copy(estimates.toDoubleArray(), 1, result, 0, estimates.size() + 1);
        return result;
        // return new double[]{validationScore / (double) n};
        //     return accumulator.createOutputStatistics(metrics);
    }

    public Collection<? extends String> getMetricNames() {
        return accumulator.metricNames();
    }
}

