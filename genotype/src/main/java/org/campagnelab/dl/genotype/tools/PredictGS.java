package org.campagnelab.dl.genotype.tools;

import edu.cornell.med.icb.util.VersionUtils;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.tools.Predict;
import org.campagnelab.dl.framework.tools.PredictArguments;
import org.campagnelab.dl.genotype.performance.BEDHelper;
import org.campagnelab.dl.genotype.performance.StatsAccumulator;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.predictions.FormatIndelVCF;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Predict tools on SSI.
 *
 * @author manuele
 */
public class PredictGS extends Predict<SegmentInformationRecords.SegmentInformation> {

    private BEDHelper bedHelper;
    private PrintWriter vcfWriter;
    private PrintWriter vcfIndelsWriter;

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


    @Override
    public PredictArguments createArguments() {
        return new PredictGSArguments();
    }

    public PredictGSArguments args() {
        return (PredictGSArguments) arguments;
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
        if (Objects.nonNull(this.vcfIndelsWriter)) this.vcfIndelsWriter.close();
        if (Objects.nonNull(this.vcfWriter)) this.vcfWriter.close();
        if (Objects.nonNull(this.bedHelper)) this.bedHelper.close();
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
        this.processAggregatedPrediction(resutsWriter, record, fullPred);
    }

    /**
     * Processes the aggregated segment prediction.
     *
     * @param resutsWriter
     * @param record
     * @param fullPred
     */
    protected void processAggregatedPrediction(PrintWriter resutsWriter,
                                               SegmentInformationRecords.SegmentInformation record, SegmentPrediction fullPred) {
        assert fullPred != null : "fullPred must not be null";
        fullPred.inspectRecord(record);

        int numOfbases = fullPred.getGenotypes().numBases();
        for (int b = 0; b < numOfbases; ) {
            // one line for each base
            ObjectList<SegmentInformationRecords.Base> basesAt = this.getBasesAt(record, b);
            int currentLocation = basesAt.get(0).getLocation();
            //TODO: to debug for the ALT empty field
            if (basesAt.get(0).getLocation() == 104364295 || basesAt.get(0).getLocation() == 104364296) {
                System.out.println("Stop here");
            }

            boolean isIndel = fullPred.isIndelPosition(basesAt.get(0));
            FormatIndelVCF format = null;
            if (isIndel) {
                Set<String> allPredicted = fullPred.predictedAllelesForIndels(b, basesAt.size());
                format = new FormatIndelVCF(fullPred.referenceAllelesForIndels(basesAt), allPredicted,
                        basesAt.get(0).getReferenceAllele().charAt(0));
            } else {
                format = new FormatIndelVCF(basesAt.get(0).getReferenceAllele(), fullPred.predictedAllelesForSnps(b),
                        basesAt.get(0).getReferenceAllele().charAt(0));
            }

            if (Objects.isNull(format))
                continue; //should never happen, but in case..

            //make an alt-allele-only set for coding
            SortedSet<String> sortedAltSet = new ObjectAVLTreeSet<String>(format.toVCF);
            sortedAltSet.remove(format.fromVCF);

            //generate alt column from alt set
            final Optional<String> optional = sortedAltSet.stream().reduce((s, s2) -> s + "," + s2);
            String altField = optional.isPresent() && !optional.get().isEmpty()? optional.get() : ".";

            //generate to column (format) from formatted predicted set
            final Optional<String> toColumnOpt = format.toVCF.stream().reduce((s, s2) -> s + "/" + s2);
            String toColumn = toColumnOpt.isPresent() ? toColumnOpt.get() : "./.";

            //get max allele length for bed file
            int maxLength = format.toVCF.stream().map(a -> a.length()).max(Integer::compareTo).orElse(0);
            maxLength = Math.max(maxLength, format.fromVCF.length());
            //TODO: if the reference is null?
            if (format.fromVCF.isEmpty() || altField.isEmpty()) {
                b += basesAt.size();
                continue;
            }
            if (args().splitIndels && isIndel) {
                // line fields: "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n";
                vcfIndelsWriter.printf(VCF_LINE, //"%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT:MC:P\t%s:%s:%f\n";
                        record.getStartPosition().getReferenceId(), //Chromosome
                        currentLocation + 1, // VCFs are 1-based
                        //ID
                        format.fromVCF, //ref
                        altField, //ALT
                        //QUAL
                        //FILTER
                        //INFO
                        PredictG.codeGT(format.toVCF, format.fromVCF, sortedAltSet),
                        toColumn,
                        fullPred.getGenotypes().probabilities[b]
                );
            } else {
                // line fields: "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n";
                vcfWriter.printf(VCF_LINE, //"%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT:MC:P\t%s:%s:%f\n";
                        record.getStartPosition().getReferenceId(), //Chromosome
                        currentLocation + 1, // VCFs are 1-based
                        //ID
                        format.fromVCF, //ref
                        altField, //ALT
                        //QUAL
                        //FILTER
                        //INFO
                        PredictG.codeGT(format.toVCF, format.fromVCF, sortedAltSet),
                        toColumn,
                        fullPred.getGenotypes().probabilities[b]
                );
                b += basesAt.size();
            }

            bedHelper.add(record.getStartPosition().getReferenceId(), currentLocation, currentLocation + maxLength, fullPred.index,
                    stats);

        }
    }

    /**
     * Returns all the bases at the same location of the base in the given position;
     *
     * @param record
     * @param position
     * @return
     */
    private ObjectList<SegmentInformationRecords.Base> getBasesAt(SegmentInformationRecords.SegmentInformation record, int position) {
        ObjectList<SegmentInformationRecords.Base> bases = new ObjectArrayList<>();
        //int b = record.getStartPosition().getLocation();
        int counter = 0;
        int detectedBaseLocation = 0;
        for (SegmentInformationRecords.Sample sample : record.getSampleList()) {
            for (SegmentInformationRecords.Base base : sample.getBaseList()) {
                if (counter++ == position || detectedBaseLocation == base.getLocation()) {
                    bases.add(base);
                    detectedBaseLocation = base.getLocation();
                }
            }
        }
        return bases;
    }

    /**
     * This method is called when we need to write the header to the results.
     *
     * @param resultsWriter
     */
    @Override
    protected void writeHeader(PrintWriter resultsWriter) {
        final String bedBasename = String.format("%s-%s-%s", modelTime, modelPrefix, testSetBasename);

        try {
            bedHelper = new BEDHelper(bedBasename);
        } catch (IOException e) {
            throw new RuntimeException("Unable to create bed file(s) to record observed regions.", e);
        }

        if (args().splitIndels) {

            final String vcfSnpsFilename = String.format("%s-%s-%s-segments-snps.vcf", modelTime, modelPrefix, testSetBasename);
            final String vcfIndelsFilename = String.format("%s-%s-%s-segments-indels.vcf", modelTime, modelPrefix, testSetBasename);

            try {
                vcfWriter = new PrintWriter(new FileWriter(vcfSnpsFilename));
                vcfWriter.append(String.format(VCF_HEADER,
                        VersionUtils.getImplementationVersion(PredictGS.class),
                        args().modelPath, args().modelName, FilenameUtils.getBaseName(args().testSet)));
            } catch (IOException e) {
                throw new RuntimeException("Unable to create VCF output file.", e);
            }

            try {
                vcfIndelsWriter = new PrintWriter(new FileWriter(vcfIndelsFilename));
                vcfIndelsWriter.append(String.format(VCF_HEADER,
                        VersionUtils.getImplementationVersion(PredictGS.class),
                        args().modelPath, args().modelName, FilenameUtils.getBaseName(args().testSet)));
            } catch (IOException e) {
                throw new RuntimeException("Unable to create VCF output file.", e);
            }
            System.out.printf("Writing VCFs and BED files: \n%s\n%s\n%s%n", vcfSnpsFilename,
                    vcfIndelsFilename, bedBasename + "-observed-regions.bed");


        } else {
            final String vcfFilename = String.format("%s-%s-%s-segments.vcf", modelTime, modelPrefix, testSetBasename);
            try {
                vcfWriter = new PrintWriter(new FileWriter(vcfFilename));
                vcfWriter.append(String.format(VCF_HEADER,
                        VersionUtils.getImplementationVersion(PredictGS.class),
                        args().modelPath, args().modelName, FilenameUtils.getBaseName(args().testSet)));
            } catch (IOException e) {
                throw new RuntimeException("Unable to create VCF output file.", e);
            }
            System.out.printf("Writing VCF and BED files: \n%s\n%s%n", vcfFilename,
                    bedBasename + "-observed-regions.bed");
        }
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
