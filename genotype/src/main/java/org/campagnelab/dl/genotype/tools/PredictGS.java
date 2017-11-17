package org.campagnelab.dl.genotype.tools;

import edu.cornell.med.icb.util.VersionUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.performance.AreaUnderTheROCCurve;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.performance.BEDHelper;
import org.campagnelab.dl.genotype.performance.StatsAccumulator;
import org.campagnelab.dl.genotype.predictions.SegmentGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.predictions.FormatIndelVCF;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * Created by mas2182 on 11/14/17.
 */
public class PredictGS extends Predict<SegmentInformationRecords.SegmentInformation> {

    private BEDHelper bedHelper;
    private PrintWriter vcfWriter;
    private static final String VCF_HEADER = "##fileformat=VCFv4.1\n" +
            "##VariationAnalysis=%s\n" +
            "##modelPath=%s\n" +
            "##modelPrefix=%s\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            "##FORMAT=<ID=MC,Number=1,Type=String,Description=\"Model Calls.\">\n" +
            "##FORMAT=<ID=P,Number=1,Type=Float,Description=\"Model proability.\">\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n";

    private static final String VCF_LINE = "%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT:MC:P\t%s:%s:%f\n";
    private String[] orderStats;
    protected StatsAccumulator stats;

    public static void main(String[] args) {

        Predict predictGS = new PredictGS();
        predictGS.parseArguments(args, "PredictGS", predictGS.createArguments());
        predictGS.execute();
    }

    public PredictGS() {
        stats = new StatsAccumulator();
        stats.initializeStats();
        orderStats = stats.createOutputHeader();
    }
    /**
     * This method is called after the test set has been observed and statistics evaluated via processPredictions.
     * It sets statistics on the whole test set, which are then written tab-delimited to a file.
     *
     * @return array of statistics to set; should be in same order as outputHeader
     */
    @Override
    protected double[] createOutputStatistics() {
        return new double[0];
    }

    /**
     * This method is called only if an output file for statistics has not been created yet. It sets the fields in
     * the header file, to then be written to the output file.
     *
     * @return array of metric names to set
     */
    @Override
    protected String[] createOutputHeader() {
        return new String[0];
    }

    /**
     * This method is called after the test set has been observed and statistics evaluated. It prints statistics
     * to the console for the operators to read.
     *
     * @param prefix The model prefix/label (e.g., bestAUC).
     */
    @Override
    protected void reportStatistics(String prefix) {

    }

    /**
     * This method is called for each record of the test set that has been predicted.
     *
     * @param resutsWriter   Where predictions can be written in tab delited format.
     * @param record
     * @param predictionList The list of predictions made for a test record. Typically one prediction for each model output.
     */
    @Override
    protected void processPredictions(PrintWriter resutsWriter, SegmentInformationRecords.SegmentInformation record, List<Prediction> predictionList) {
        SegmentPrediction fullPred = (SegmentPrediction) domainDescriptor.aggregatePredictions(record, predictionList);
        assert fullPred != null : "fullPref must not be null";

        System.out.println(fullPred.getGenotypes());

        //fullPred.inspectRecord(record);
        int bases = fullPred.getGenotypes().numBases();
        int startPosition = fullPred.getStartPosition();
        for (int b =0; b < fullPred.getGenotypes().numBases(); b++) {
            // one line for each base
            //FormatIndelVCF format = new FormatIndelVCF(fullPred.getReferenceId(), fullPred.(), fullPred.predictedFrom.charAt(0));


              /*vcfWriter.printf(VCF_LINE, fullPred.getReferenceId(), fullPred.getStartPosition(), fullPred.getEndPosition(),
                    format.fromVCF, altField, codeGT(format.toVCF, format.fromVCF, sortedAltSet), toColumn,
                    fullPred.overallProbability);
              */
        }

    }

    /**
     * This method is called when we need to write the header to the results.
     *
     * @param resultsWriter
     */
    @Override
    protected void writeHeader(PrintWriter resultsWriter) {
        final String vcfFilename = String.format("%s-%s-%s-genotypes.vcf", modelTime, modelPrefix, testSetBasename);
        final String bedBasename = String.format("%s-%s-%s", modelTime, modelPrefix, testSetBasename);

        try {
            vcfWriter = new PrintWriter(new FileWriter(vcfFilename));
        } catch (IOException e) {
            throw new RuntimeException("Unable to create VCF output file.", e);
        }

        vcfWriter.append(String.format(VCF_HEADER,
                VersionUtils.getImplementationVersion(PredictG.class),
                args().modelPath, args().modelName, FilenameUtils.getBaseName(args().testSet)));
        try {
            bedHelper = new BEDHelper(bedBasename);
        } catch (IOException e) {
            throw new RuntimeException("Unable to create bed file(s) to record observed regions.", e);
        }
        System.out.printf("Writing VCF and BED files: \n%s\n%s%n", vcfFilename, bedBasename + "-observed-regions.bed");

    }

    /**
     * This method must allocate any statistic calculator necessary to evaluate performance on the test set and store
     * these instances in fields of the sub-class.
     *
     * @param prefix The model prefix/label (e.g., bestAUC).
     */
    @Override
    protected void initializeStats(String prefix) {
        stats = new StatsAccumulator();
        stats.initializeStats();
        orderStats = stats.createOutputHeader();
    }
}
