package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.performance.StatsAccumulator;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.PrintWriter;
import java.util.Arrays;
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
    /**
     * We estimate the AUC of a correct prediction on a variant with respect to everything else (e.g., incorrect
     * on variant or reference).
     */
    private AreaUnderTheROCCurve aucLossCalculator;
    private double auc;
    private double[] confidenceInterval95;

    @Override
    public PredictArguments createArguments() {
        return new PredictGArguments();
    }

    public static void main(String[] args) {

        Predict predict = new PredictG();
        predict.parseArguments(args, "PredictG", predict.createArguments());
        predict.execute();
    }


    protected StatsAccumulator stats;

    @Override
    protected void writeHeader(PrintWriter resutsWriter) {
        resutsWriter.append("index\tpredictionCorrect01\ttrueGenotypeCall\tpredictedGenotypeCall\tprobabilityIsCalled\tcorrectness\tregion\tisVariant").append("\n");
    }

    @Override
    protected void initializeStats(String prefix) {
        stats = new StatsAccumulator();
        stats.setNumVariantsExpected(args().numVariantsExpected);
        stats.initializeStats();
        aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
    }

    String[] orderStats = {"Accuracy", "Recall", "Precision",
            "F1", "NumVariants", "Concordance"};

    @Override
    protected double[] createOutputStatistics() {
        DoubleList values = new DoubleArrayList();

        values.addAll(DoubleArrayList.wrap(stats.createOutputStatistics(orderStats)));
        auc = aucLossCalculator.evaluateStatistic();
        confidenceInterval95 = aucLossCalculator.confidenceInterval95();

        values.add(auc);
        values.add(confidenceInterval95[0]);
        values.add(confidenceInterval95[1]);

        return values.toDoubleArray();
    }

    @Override
    protected String[] createOutputHeader() {


        ObjectArrayList<String> values = new ObjectArrayList();

        values.addAll(it.unimi.dsi.fastutil.objects.ObjectArrayList.wrap(orderStats));
        values.add("AUC");
        values.add("[AUC95");
        values.add("AUC95]");
        return values.toArray(new String[0]);
    }

    @Override
    protected void reportStatistics(String prefix) {
        stats.reportStatistics(prefix);
        System.out.printf("AUC = %f [%f-%f]%n", auc,
                confidenceInterval95[0], confidenceInterval95[1]);
        System.out.println("Printable: " + Arrays.toString(createOutputStatistics()));
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
        if (GenotypeHelper.isNoCall(fullPred.predictedGenotype)) {
            System.out.printf("preventing no call from being interpreted as a variant: %s %s %n", fullPred.predictedGenotype, record.getReferenceBase());
            fullPred.isVariant = false;
        }
        boolean correct = fullPred.isCorrect();
        //remove dangling commas
        String correctness = correct ? "correct" : "wrong";

        if (filterHet(args(), fullPred) &&
                filterVariant(args(), fullPred) &&
                doOuptut(correctness, args(), fullPred.overallProbability)) {
            resultWriter.printf("%d\t%d\t%s\t%s\t%f\t%s\t%s:%s\t%s\n",
                    fullPred.index, (correct ? 1 : 0),
                    fullPred.trueGenotype, fullPred.predictedGenotype,
                    fullPred.isVariantProbability, correctness, record.getReferenceId(), record.getPosition() + 1,
                    fullPred.isVariant ? "variant" : "-");
            if (args().filterMetricObservations) {
                stats.observe(fullPred);
                observeForAUC(fullPred);
            }
        }
        if (!args().filterMetricObservations) {
            stats.observe(fullPred);
            observeForAUC(fullPred);
        }

    }

    private void observeForAUC(GenotypePrediction fullPred) {
        if (fullPred.isVariant()) {
            aucLossCalculator.observe(fullPred.overallProbability,  fullPred.isCorrect() ? 1 : -1);
        }
    }

    private boolean filterVariant(PredictGArguments args, GenotypePrediction fullPred) {
        if (args().onlyVariants) {
            return fullPred.isVariant();
        } else {
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
