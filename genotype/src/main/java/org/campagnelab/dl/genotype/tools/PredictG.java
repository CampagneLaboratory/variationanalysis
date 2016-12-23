package org.campagnelab.dl.genotype.tools;


import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.performance.StatsAccumulator;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Remi Torracinta
 *         Created by rct66 on 12/7/16.
 */
public class PredictG extends Predict<BaseInformationRecords.BaseInformation> {

    private AreaUnderTheROCCurve aucLossCalculator;

    @Override
    public PredictArguments createArguments() {
        return new PredictGArguments();
    }

    public static void main(String[] args) {

        Predict predict = new PredictG();
        predict.parseArguments(args, "PredictG", predict.createArguments());
        predict.execute();
    }


    protected StatsAccumulator stats = new StatsAccumulator();

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        resutsWriter.append("index\tpredictionCorrect01\ttrueGenotypeCall\tpredictedGenotypeCall\tprobabilityIsCalled\tcorrectness").append("\n");
    }

    @Override
    protected void initializeStats(String prefix) {
        stats.initializeStats();
        aucLossCalculator = new AreaUnderTheROCCurve();
    }


    @Override
    protected double[] createOutputStatistics() {
        return stats.createOutputStatistics("Accuracy", "Recall", "Precision",
                "F1", "NumVariants", "Concordance", "score");
    }

    @Override
    protected String[] createOutputHeader() {
        return stats.createOutputHeader();
    }

    @Override
    protected void reportStatistics(String prefix) {
        stats.reportStatistics(prefix);
    }

    public PredictGArguments args() {
        return (PredictGArguments) arguments;
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, BaseInformationRecords.BaseInformation record, List<Prediction> predictionList) {


        GenotypePrediction fullPred = (GenotypePrediction) domainDescriptor.aggregatePredictions(predictionList);
        fullPred.inspectRecord(record);
        long trueAlleleLength = fullPred.trueAlleles().stream().map(String::length).distinct().count();
        if (!args().scoreIndels && (fullPred.isIndel || trueAlleleLength > 1)) {
            // reduce A---A/ATTTA to A/A
            String trimmedGenotype = fullPred.trueAlleles().stream().map(s -> Character.toString(s.charAt(0))).collect(Collectors.joining("/"));
            fullPred.trueGenotype = trimmedGenotype;
            // no longer an indel, and now matching reference:
            fullPred.isIndel = false;
            fullPred.isVariant = false;
        }
        boolean correct = fullPred.isCorrect();
        //remove dangling commas
        String correctness = correct ? "correct" : "wrong";

        if (filterHet(args(), fullPred) &&
                filterVariant(args(), fullPred) &&
                doOuptut(correctness, args(), fullPred.overallProbability)) {
            resultWriter.printf("%d\t%d\t%s\t%s\t%f\t%s\n",
                    fullPred.index, (correct ? 1 : 0), fullPred.trueGenotype, fullPred.predictedGenotype, fullPred.overallProbability, correctness);
            if (args().filterMetricObservations) {
                stats.observe(fullPred);
            }
        }
        if (!args().filterMetricObservations) {
            stats.observe(fullPred);
        }

    }

    private boolean filterVariant(PredictGArguments args, GenotypePrediction fullPred) {
        if (args().onlyVariants) {
            return fullPred.isVariant();
        }
        else {
            return true;
        }
    }

    boolean filterHet(PredictGArguments args, GenotypePrediction fullPred) {
        Set<String> alleles = fullPred.predictedAlleles();
        switch (args.showFilter) {
            case HET:
                return (alleles.size() == 2);
            case HOM:
                return alleles.size() == 1;
            default:
                return true;
        }
    }

    /**
     * Apply filters and decide if a prediction should be written to the output.
     *
     * @param correctness
     * @param args
     * @param pMax
     * @return
     */
    protected boolean doOuptut(String correctness, PredictArguments args, double pMax) {
        if (args.correctnessFilter != null) {
            if (!correctness.equals(args.correctnessFilter)) {
                return false;
            }
        }
        if (pMax < args().pFilterMinimum || pMax > args().pFilterMaximum) {
            return false;
        }
        return true;
    }


}
