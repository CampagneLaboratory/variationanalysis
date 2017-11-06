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
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Generate SBI files starting from a variant map.
 */
public class SBISimulator extends AbstractTool<SBISimulatorArguments> {

    static private Logger LOG = LoggerFactory.getLogger(SBISimulator.class);
    private Random r = new Random();

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
        try (RecordWriter writer = new RecordWriter(args().outputFilename)) {
            final VariantMapHelper helper = new VariantMapHelper(args().inputFile);
            List<String> chromosomes = this.chromosomesForSBI(helper);
            chromosomes.parallelStream().limit(args().readN).forEach((String chromosome) -> {
                System.out.println("Chrom: " + chromosome);
                for (ObjectIterator<Variant> it = helper.getAllVariants(chromosome); it.hasNext(); ) {
                    Variant variant = it.next();
                    for (Variant.FromTo trueAllele : variant.trueAlleles) {
                        String[] counts;
                        String genotype = String.format("%s/%s", trueAllele.getFrom(), trueAllele.getTo());
                        if (trueAllele.getFrom().equals(trueAllele.getTo())) {
                            counts = new String[]{
                                    String.format(countFormat, trueAllele.getFrom(), trueAllele.getTo(), generateCounts(), generateCounts())};
                        } else {
                            counts = new String[]{
                                    String.format(countFormat, trueAllele.getFrom(), trueAllele.getFrom(), generateCounts(), generateCounts()),
                                    String.format(countFormat, trueAllele.getFrom(), trueAllele.getTo(), generateCounts(), generateCounts())};
                        }
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

}
