package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
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
    private char[] bases = new char[]{'A', 'C', 'T', 'G', 'N'};
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
        try {
            final RecordWriter writer = new RecordWriter(args().outputFilename);
            final VariantMapHelper helper = new VariantMapHelper(args().inputFile);
            List<String> chromosomes = this.chromosomesForSBI(helper);
            chromosomes.parallelStream().limit(args().readN).forEach((String chromosome) -> {
                if (args().verbose) System.out.println("Chrom: " + chromosome);
                for (ObjectIterator<Variant> it = helper.getAllVariants(chromosome); it.hasNext(); ) {
                    Variant variant = it.next();
                    String trueGenotype = GenotypeHelper.fromFromTos(variant.trueAlleles);
                    Set<String> counts = new HashSet<>();
                    Set<String> keys = new HashSet<>();
                    if (GenotypeHelper.getAlleles(trueGenotype).size() > 1) {
                        String keepFrom = variant.referenceBase;
                        for (Variant.FromTo trueAllele : variant.trueAlleles) {                         
                            if (!trueAllele.getFrom().equals(variant.referenceBase))
                                keepFrom = trueAllele.getFrom(); //use the one that does not match the reference, if any
                        }
                        for (Variant.FromTo trueAllele : variant.trueAlleles) {
                            calculateCounts(keepFrom, (trueAllele.getTo().length() != keepFrom.length()) ? keepFrom : trueAllele.getTo(),
                                    counts, keys, trueGenotype);
                        }
                    }  else {
                        for (Variant.FromTo trueAllele : variant.trueAlleles) {
                            calculateCounts(trueAllele.getFrom(), trueAllele.getTo(), counts, keys, trueGenotype);
                        }
                    }
                    fillUpCounts(counts, keys, variant.referenceBase);
                    if (args().verbose)
                        System.out.println("Genotype=" + trueGenotype + ", counts: " +
                                Arrays.toString(counts.toArray(new String[counts.size()])));
                    BaseInformationRecords.BaseInformation record = makeRecord(variant.referenceIndex, chromosome,
                            variant.position, trueGenotype, counts.toArray(new String[counts.size()]));
                    addRecord(writer, record);
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
    private BaseInformationRecords.BaseInformation makeRecord(int refIndex, String refId, int position, String genotype, String... countCreations) {
        BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        builder.setTrueGenotype(genotype);
        builder.setReferenceIndex(refIndex);
        builder.setPosition(position);
        builder.setReferenceId(refId);
        BaseInformationRecords.SampleInfo.Builder sample = BaseInformationRecords.SampleInfo.newBuilder();
        String referenceBase = "N";
        for (String countCreationInstruction : countCreations) {
            BaseInformationRecords.CountInfo.Builder countBuilder = BaseInformationRecords.CountInfo.newBuilder();
            String tokens[] = countCreationInstruction.split("[/=+]");
            assert tokens.length == 4 :
                    "count creation instruction must have four arguments: ref/to=forward+reverse, was " + countCreationInstruction;
            final String from = tokens[0];
            countBuilder.setFromSequence(from);
            referenceBase = Character.toString(from.charAt(0));
            builder.setReferenceBase(referenceBase);
            final String token = tokens[1];
            countBuilder.setIsCalled(true);
            countBuilder.setToSequence(token);
            countBuilder.setMatchesReference(from.equals(token));
            countBuilder.setGenotypeCountForwardStrand(Integer.parseInt(tokens[2]));
            countBuilder.setGenotypeCountReverseStrand(Integer.parseInt(tokens[3]));
            populateWithFrequencies(countBuilder);
            if (from.length() > 0) {
                countBuilder.setIsIndel(true);
            }
            sample.addCounts(countBuilder);
        }
        builder.setGenomicSequenceContext(referenceBase);
        sample.setFormattedCounts(FormatterCountHelper.format(sample));
        builder.addSamples(sample);
        return builder.build();
    }

    private void populateWithFrequencies(BaseInformationRecords.CountInfo.Builder countBuilder) {
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

    private void calculateCounts(String from, String to, Set<String> allCounts, Set<String> keys, String trueGenotype) {
       /* String[] inGenotype = trueGenotype.split("\\|", 2);
                if ((from.equals(inGenotype[0]) || to.equals(inGenotype[1]))
                        && !keys.contains(from + "/" + to)) {
                    allCounts.add(String.format(countFormat, from, to, generateCounts(), generateCounts()));
                    keys.add(from + "/" + to);
                }  */
        allCounts.add(String.format(countFormat, from, to, generateCounts(), generateCounts()));
        keys.add(from + "/" + to);
    }

    private void fillUpCounts(Set<String> allCounts, Set<String> keys, String referenceBase) {
        for (char base : bases) {
            String key = referenceBase + "/" + base;
            if (!keys.contains(key)) {
                allCounts.add(String.format(countFormat, referenceBase, base, 0, 0));
                keys.add(key);
            }
        }
    }
}
