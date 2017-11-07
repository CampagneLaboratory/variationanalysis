package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.segments.FormatterCountHelper;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.util.Variant;
import org.campagnelab.goby.util.VariantMapHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.*;

/**
 * Generate SBI files starting from a variant map.
 */
public class SBISimulator extends AbstractTool<SBISimulatorArguments> {

    static private Logger LOG = LoggerFactory.getLogger(SBISimulator.class);
    private Random r = new Random();
    private char[] bases = new char[]{'A','C','T','G','N'};
    private String countFormat = "%s/%s=%d+%d";

    public static void main(String[] args) {
        SBISimulator tool = new SBISimulator();
        tool.parseArguments(args, "SBISimulator", tool.createArguments());
        tool.execute();
    }

    @Override
    public SBISimulatorArguments createArguments() {
        return new SBISimulatorArguments();
    }

    @Override
    public void execute() {
        String countFormat = "%s/%s=%d+%d";
        try {
            final RecordWriter writer = new RecordWriter(args().outputFilename);
            final VariantMapHelper helper = new VariantMapHelper(args().inputFile);
            List<String> chromosomes = this.chromosomesForSBI(helper);
            chromosomes.parallelStream().limit(args().readN).forEach((String chromosome) -> {
                System.out.println("Chrom: " + chromosome);
                for (ObjectIterator<Variant> it = helper.getAllVariants(chromosome); it.hasNext(); ) {
                    Variant variant = it.next();
                    for (Variant.FromTo trueAllele : variant.trueAlleles) {
                        String[] counts = fillUpCounts(trueAllele, variant.referenceBase);
                        String genotype = String.format("%s/%s", trueAllele.getFrom(), trueAllele.getTo());
                        System.out.println("Genotype=" + genotype + ", counts: " + Arrays.toString(counts));
                        BaseInformationRecords.BaseInformation record = makeRecord(variant.referenceIndex, variant.position, genotype, counts);
                        addRecord(writer, record);
                    }
                }
            });
            writer.close();
        } catch (IOException e) {
            LOG.error("Unable to locate the variant map.");
            throw new IllegalArgumentException("Unable to locate the variant map.");
        } catch (ClassNotFoundException e) {
            LOG.error("Unable to load the variant map.");
            throw new IllegalStateException("Unable to load the variant map.");
        }
    }

    /**
     * Detect the chromosomes to use for the output SBI
     *
     * @param helper
     * @return
     */
    private List<String> chromosomesForSBI(final VariantMapHelper helper) {
        List<String> chroms = null;
        if (args().chromosome != null) {
            //need to override what is in the map
            chroms = new ObjectArrayList<String>(1);
            chroms.add(args().chromosome);
        } else {
            chroms = new ObjectArrayList<String>(helper.size());
            ObjectIterator<String> it = helper.getAllChromosomes();
            while (it.hasNext()) chroms.add(it.next());
        }
        return chroms;
    }

    // format of count creation instruction is from/to=10+12
    private BaseInformationRecords.BaseInformation makeRecord(int refIndex, int position, String genotype, String... countCreations) {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setTrueGenotype(genotype);
        builder.setReferenceIndex(refIndex);
        builder.setPosition(position);

        BaseInformationRecords.SampleInfo.Builder sample = BaseInformationRecords.SampleInfo.newBuilder();
        for (String countCreationInstruction : countCreations) {
            BaseInformationRecords.CountInfo.Builder countBuilder = BaseInformationRecords.CountInfo.newBuilder();
            String tokens[] = countCreationInstruction.split("[/=+]");
            assert tokens.length == 4 :
                    "count creation instruction must have four arguments: ref/to=forward+reverse, was " + countCreationInstruction;
            final String from = tokens[0];
            countBuilder.setFromSequence(from);
            builder.setReferenceBase(Character.toString(from.charAt(0)));
            final String token = tokens[1];
            countBuilder.setToSequence(token);
            countBuilder.setMatchesReference(from.equals(token));
            countBuilder.setGenotypeCountForwardStrand(Integer.parseInt(tokens[2]));
            countBuilder.setGenotypeCountReverseStrand(Integer.parseInt(tokens[3]));
            BaseInformationRecords.NumberWithFrequency.Builder builderN = BaseInformationRecords.NumberWithFrequency.newBuilder();
            builderN.setFrequency(1);
            builderN.setNumber(1);
            BaseInformationRecords.NumberWithFrequency frequency = builderN.build();
            countBuilder.addQueryPositions(frequency);
            countBuilder.addNumVariationsInReads(frequency);
            countBuilder.addDistancesToReadVariationsForwardStrand(frequency);
            countBuilder.addDistancesToReadVariationsReverseStrand(frequency);
            countBuilder.addDistanceToStartOfRead(frequency);
            countBuilder.addDistanceToEndOfRead(frequency);
            countBuilder.addReadMappingQualityForwardStrand(frequency);
            countBuilder.addReadIndicesForwardStrand(frequency);
            countBuilder.addReadIndicesReverseStrand(frequency);
            countBuilder.addQualityScoresForwardStrand(frequency);
            countBuilder.addQualityScoresReverseStrand(frequency);
            countBuilder.addTargetAlignedLengths(frequency);
            countBuilder.addQueryAlignedLengths(frequency);
            countBuilder.addReadMappingQualityForwardStrand(frequency);
            countBuilder.addReadMappingQualityReverseStrand(frequency);
            if (from.contains("-") || token.contains("-")) {
                countBuilder.setIsIndel(true);
            }
            sample.addCounts(countBuilder);
        }
        sample.setFormattedCounts(FormatterCountHelper.format(sample));
        builder.addSamples(sample);
        return builder.build();
    }

    private void addRecord(RecordWriter writer, BaseInformationRecords.BaseInformation record) {
        try {
            writer.writeRecord(record);
        } catch (IOException e) {
            LOG.error("Unable to write the record.", e);
        }
    }

    private int generateCounts() {
        return r.nextInt(100);
    }

    private String[] fillUpCounts(final Variant.FromTo allele, String referenceBase) {
        Set<String> allCounts = new HashSet<>();
        ObjectArrayList<String> counts = new ObjectArrayList<>();
        if (allele.getFrom().equals(allele.getTo())) {
            counts.add(String.format(countFormat, allele.getFrom(), allele.getTo(), generateCounts(), generateCounts()));
            allCounts.add(allele.getFrom()+"/"+allele.getTo());
        } else {
            counts.add(String.format(countFormat, allele.getFrom(), allele.getFrom(), generateCounts(), generateCounts()));
            counts.add(String.format(countFormat, allele.getFrom(), allele.getTo(), generateCounts(), generateCounts()));
            allCounts.add(allele.getFrom()+"/"+allele.getFrom());
            allCounts.add(allele.getFrom()+"/"+allele.getTo());
        }
        for (char base : bases){
            String key = referenceBase +"/" + base;
            if (!allCounts.contains(key)) {
                counts.add(String.format(countFormat, referenceBase, base, 0, 0));
                allCounts.add(key);
            }
        }

        return counts.toArray(new String[counts.size()]);
    }

}
