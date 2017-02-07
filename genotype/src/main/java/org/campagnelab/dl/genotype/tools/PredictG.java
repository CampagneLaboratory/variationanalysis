package org.campagnelab.dl.genotype.tools;

import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.performance.StatsAccumulator;
import org.campagnelab.dl.genotype.predictions.FormatIndelVCF;
import org.campagnelab.dl.genotype.predictions.GenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.alignments.ConcatSortedAlignmentReader;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Example of Predict implementation. This class performs predictions with a model trained by TrainModelS.
 *
 * @author Remi Torracinta
 * @author Fabien Campagne
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
    private PrintWriter bedWriter;
    private PrintWriter vcfWriter;


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
        final String vcfFilename = String.format("%s-%s-%s-genotypes.vcf", modelTime, modelPrefix, testSetBasename);
        final String bedFilename = String.format("%s-%s-%s-observed-regions.bed", modelTime, modelPrefix, testSetBasename);

        try {
            vcfWriter = new PrintWriter(new FileWriter(vcfFilename));
        } catch (IOException e) {
            throw new RuntimeException("Unable to create VCF output file.", e);
        }

        if (args().outputFormat == PredictGArguments.OutputFormat.VCF) {
            vcfWriter.append(String.format(VCF_HEADER,
                    VersionUtils.getImplementationVersion(PredictG.class),
                    args().modelPath, args().modelName));
        } else {
            resutsWriter.append("index\tpredictionCorrect01\ttrueGenotypeCall\tpredictedGenotypeCall\tprobabilityIsCalled\tcorrectness\tregion\tisVariant").append("\n");
        }
        try {
            bedWriter = new PrintWriter(new FileWriter(bedFilename));
        } catch (IOException e) {
            throw new RuntimeException("Unable to create bed file to record observed regions.", e);
        }

        System.out.printf("Writing VCF and BED files: \n%s\n%s%n", vcfFilename, bedFilename);

    }

    @Override
    protected void initializeStats(String prefix) {
        stats = new StatsAccumulator();
        stats.setNumVariantsExpected(args().numVariantsExpected);
        stats.initializeStats();
        aucLossCalculator = new AreaUnderTheROCCurve(args().numRecordsForAUC);
    }

    String[] orderStats = {"Accuracy", "Recall", "Precision", "F1", "NumVariants", "Concordance",
            "Accuracy_Indels", "Recall_Indels", "Precision_Indels", "F1_Indels",
            "Accuracy_SNPs", "Recall_SNPs", "Precision_SNPs", "F1_SNPs",
    };

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

    private static final String VCF_HEADER = "##fileformat=VCFv4.1\n" +
            "##VariationAnalysis=%s\n" +
            "##modelPath=%s\n" +
            "##modelPrefix=%s\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            "##FORMAT=<ID=MC,Number=1,Type=String,Description=\"Model Calls.\">\n" +
            "##FORMAT=<ID=P,Number=1,Type=Float,Description=\"Model proability.\">\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n";

    private static final String VCF_LINE = "%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT:MC:P\t%s:%s:%f\n";

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
        bedWriter.close();
        vcfWriter.close();
    }

    public PredictGArguments args() {
        return (PredictGArguments) arguments;
    }

    @Override
    protected void processPredictions(PrintWriter resultWriter, BaseInformationRecords.BaseInformation record, List<Prediction> predictionList) {

        if (coverage(record) < args().minimumCoverage) {
            return;
        }
        GenotypePrediction fullPred = (GenotypePrediction) domainDescriptor.aggregatePredictions(predictionList);
        fullPred.inspectRecord(record);
        long trueAlleleLength = fullPred.trueAlleles().stream().map(String::length).distinct().count();
        if (!args().scoreIndels && (fullPred.isIndel || trueAlleleLength > 1)) {
            // reduce A---A/ATTTA to A/A
            String trimmedGenotype = fullPred.trueAlleles().stream().map(s -> Character.toString(s.charAt(0))).collect(Collectors.joining("/"));
            fullPred.trueGenotype = trimmedGenotype;
            fullPred.trueFrom = record.getReferenceBase();
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
        // obtain isVariant from the gold-standard, not from the prediction.
        boolean isVariant = record.getSamples(0).getIsVariant();
        final boolean isPredictedVariant = GenotypeHelper.isVariant(fullPred.predictedAlleles(), record.getReferenceBase());

        if (filterHet(args(), fullPred) &&
                filterVariant(args(), fullPred) &&
                doOuptut(correctness, args(), fullPred.overallProbability)) {
            switch (args().outputFormat) {
                case TABLE:
                    resultWriter.printf("%d\t%d\t%s\t%s\t%f\t%s\t%s:%s\t%s\t\n",
                            fullPred.index, (correct ? 1 : 0),
                            fullPred.trueGenotype, fullPred.predictedGenotype,
                            fullPred.isVariantProbability, correctness, record.getReferenceId(), record.getPosition() + 1,
                            isVariant ? "variant" : "-");
                    break;
                case VCF:
                    //TODO: improve logic so that a heterozygous SNP/Indel at the same position is handled. currently, the snp isn't adjusted to correspond to the indel's from field.
                    //generated vcf formatted indel
                    FormatIndelVCF format = new FormatIndelVCF(fullPred.predictedFrom,fullPred.predictedAlleles(),fullPred.predictedFrom.charAt(0));

                    //get max allele length for bed file
                    int maxLength = format.toVCF.stream().map(a -> a.length()).max(Integer::compareTo).orElse(0);
                    maxLength = Math.max(maxLength,format.fromVCF.length());

                    //make an alt-allele-only set for coding
                    SortedSet<String> sortedAltSet = new ObjectAVLTreeSet<String>(format.toVCF);
                    sortedAltSet.remove(format.fromVCF);

                    //generate alt column from alt set
                    final Optional<String> optional = sortedAltSet.stream().reduce((s, s2) -> s + "," + s2);
                    String altField = optional.isPresent() ? optional.get() : ".";

                    //generate to column (format) from formatted predicted set
                    final Optional<String> toColumnOpt = format.toVCF.stream().reduce((s, s2) -> s + "/" + s2);
                    String toColumn = toColumnOpt.isPresent() ? toColumnOpt.get() : "./.";


                    if (sortedAltSet.size() >= 1) {
                        // only append to VCF if there is at least one alternate allele:
                        // NB: VCF format is one-based.
                        vcfWriter.printf(VCF_LINE, record.getReferenceId(), record.getPosition() + 1,
                                format.fromVCF, altField, codeGT(format.toVCF, format.fromVCF, sortedAltSet), toColumn, fullPred.isVariantProbability);
                    }
                    // NB: bed format is zero-based.
                    bedWriter.printf("%s\t%d\t%d\t%d\n", record.getReferenceId(), record.getPosition(), record.getPosition() + maxLength, fullPred.index);

                    break;
            }
            if (args().filterMetricObservations) {

                stats.observe(fullPred, isVariant, isPredictedVariant);
                observeForAUC(fullPred, isVariant);
            }
        }
        if (!args().filterMetricObservations) {
            stats.observe(fullPred, fullPred.isVariant(), isPredictedVariant);
            observeForAUC(fullPred, isVariant);
        }

    }

    private int coverage(BaseInformationRecords.BaseInformation record) {
        int coverage = 0;
        for (BaseInformationRecords.CountInfo count : record.getSamples(0).getCountsList())
            coverage += count.getGenotypeCountForwardStrand() + count.getGenotypeCountReverseStrand();
        return coverage;
    }

    /**
     * @param to
     * @param from
     * @param altSet
     * @return
     */
    public static String codeGT(Set<String> to, String from, SortedSet<String> altSet) {
        IntArrayList codedAlleles = new IntArrayList();

        String result = "";
        for (String allele : to) {
            if (from.equals(allele)) {
                codedAlleles.add(0);
            }
            int altIndex = 1;
            for (String altAllele : altSet) {
                if (altAllele.equals(allele)) {
                    codedAlleles.add(altIndex);
                }
                altIndex += 1;
            }
        }
        if (codedAlleles.size() > 0) {
            Collections.sort(codedAlleles);
            result = "";
            boolean first = true;
            for (int code : codedAlleles) {
                if (!first) {
                    result += "/";
                }
                result += Integer.toString(code);
                first = false;
            }
            return result;
        } else {
            return "./.";
        }
    }

    private void observeForAUC(GenotypePrediction fullPred, boolean isVariant) {
        if (isVariant) {
            aucLossCalculator.observe(fullPred.overallProbability, fullPred.isCorrect() ? 1 : -1);
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
