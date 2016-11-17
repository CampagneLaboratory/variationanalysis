package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.somatic.intermediaries.*;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.intermediaries.SimulationStrategy;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.intermediaries.SimulationStrategyImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

/**
 * The mutator object iterates over a file and creates additional copies of every record, where
 * some copies are mutated, and some are not. In constrast to Mutator, Mutator2 observes the sum of counts in the
 * second sample, discards the counts in the second sample, and recreates them using a simulation protocol
 * that either simulates random sampling from the same genotype, or from a second tumor genotype derived from germline.
 * <p>
 *
 * @author Fabien Campagne
 */
public class Mutator2 extends AbstractTool<Mutator2Arguments> {
    private static final int CHUNK_SIZE = 10000;
    private static final int NUM_SIMULATED_RECORD_PER_DATUM = 2;
    static private Logger LOG = LoggerFactory.getLogger(Mutator2.class);

    final String[] STRING = new String[]{"A", "T", "C", "G"};
    Random rand;
    final boolean MUTATE = true;
    SimulationStrategy strategy;
    boolean trioUnchecked = true;
    int numCanonical = 0;
    int numRecordsTotal = 0;
    final double deltaSmall = 0.0;
    final double deltaBig = 1.0;
    final int seed = 2323;
    final int seed2 = 2348999;
    public static void main(String[] args) throws IOException {



        //new ParquetPrinter(args[0]).print();
//        if (args.length < 2) {
//            System.err.println("usage: input.sbi mutated-randomized-filename");
//            System.exit(1);
//        }
        Mutator2 m = new Mutator2(args);
        m.execute();
        //  new ParquetPrinter(args[2]).print();

        System.out.println("Fraction of non-canonical:" + ((float)1-((float)m.numCanonical/(float)m.numRecordsTotal)));
    }

    public Mutator2(String [] args) {
        setSeed(seed);
        this.parseArguments(args, "Mutator2", this.createArguments());
        strategy = new SimulationStrategyImpl(deltaSmall,deltaBig,args().heteroHeuristic,seed2,args().canonThreshold);
    }

    public void setSeed(int seed) {

        rand = new XorShift1024StarRandom(seed);
    }

    //deprectated method should not be used, arguments obatined from jcommander
//    public Mutator2(int deltaSmall, int deltaBig, int zygHeuristic) {
//        this.deltaSmall = deltaSmall;
//        this.deltaBig = deltaBig;
//        this.zygHeuristic = zygHeuristic;
//        setSeed(2323);
//    }


    public void execute() {
        try {
            RecordReader reader = new RecordReader(args().inputFile);
            RecordWriter writer = new RecordWriter(args().outputFile);

            //set up logger
            ProgressLogger pgReadWrite = new ProgressLogger(LOG);
            pgReadWrite.itemsName = "mutation";
            pgReadWrite.expectedUpdates = reader.getTotalRecords();
            pgReadWrite.displayFreeMemory = true;
            pgReadWrite.start();
            SimulationCharacteristics sim = new SimulationCharacteristics();
            int maxProcess = Integer.MAX_VALUE;
            int iteration = 0;
            for (BaseInformationRecords.BaseInformation base : reader) {
                iteration++;
                sim.observe(base);
                if (sim.size() >= CHUNK_SIZE) {
                    sim.batchIsComplete();
                    processBatch(sim, writer);
                    sim.clear();
                }

                pgReadWrite.lightUpdate();
                if (iteration > maxProcess) {
                    break;
                }
            }
            processBatch(sim, writer);

            pgReadWrite.stop();
            reader.close();
            writer.close();
        } catch (IOException e) {
            System.err.println("Unable to load or write files. Check command line arguments.");
        }
    }

    /**
     * Process one batch of records.
     *
     * @param sim
     * @param writer
     */
    private void processBatch(SimulationCharacteristics sim, RecordWriter writer) throws IOException {
        Iterator<BaseInformationRecords.BaseInformation> iterator = sim.iterator();
        ObjectArrayList<BaseInformationRecords.BaseInformation> shufflingList = new ObjectArrayList<>();
        while (iterator.hasNext()) {
            BaseInformationRecords.BaseInformation record = iterator.next();
            if (trioUnchecked && record.getSamplesCount() > 2){
                strategy = new SimulationStrategyImplTrio(deltaSmall, deltaBig, args().heteroHeuristic, seed2, args().canonThreshold);
                trioUnchecked = false;
            }
            shufflingList.add(strategy.mutate(false, record, record.getSamples(0), record.getSamples(1), sim));
            numRecordsTotal++;

            for (int i = 0; i < NUM_SIMULATED_RECORD_PER_DATUM; i++) {
                BaseInformationRecords.BaseInformation possiblyMutated = strategy.mutate(true, record, record.getSamples(0), record.getSamples(1), sim);
                if (possiblyMutated.getMutated()) {
                    shufflingList.add(possiblyMutated);
                    if (i==0){
                        numCanonical++;
                    }
                } else {
                    break;
                }
            }
        }
        Collections.shuffle(shufflingList);
        double numMutated = 0;
        for (BaseInformationRecords.BaseInformation record : shufflingList) {

            writer.writeRecord(record);
            numMutated += record.getMutated() ? 1 : 0;
        }
        System.out.printf("Ratio of mutated to total record (0-1): %f%n", numMutated / shufflingList.size());
        System.out.flush();
        shufflingList.clear();
    }

    @Override
    public Mutator2Arguments createArguments() {
        return new Mutator2Arguments();
    }

    public static String regenerateFormattedCounts(BaseInformationRecords.SampleInfoOrBuilder sample, String mutatedAllele) {
        int a = sample.getCounts(0).getGenotypeCountReverseStrand() + sample.getCounts(0).getGenotypeCountForwardStrand();
        int t = sample.getCounts(1).getGenotypeCountReverseStrand() + sample.getCounts(1).getGenotypeCountForwardStrand();
        int c = sample.getCounts(2).getGenotypeCountReverseStrand() + sample.getCounts(2).getGenotypeCountForwardStrand();
        int g = sample.getCounts(3).getGenotypeCountReverseStrand() + sample.getCounts(3).getGenotypeCountForwardStrand();
        int n = sample.getCounts(4).getGenotypeCountReverseStrand() + sample.getCounts(4).getGenotypeCountForwardStrand();
        String fb;
        try {
            fb = sample.getFormattedCounts().split(" ")[8];
        } catch (ArrayIndexOutOfBoundsException e) {
            fb = "n/a";
        }
        int numIndels = sample.getCountsCount() - 5;
        int[] indels = new int[numIndels];
        for (int i = 5; i < numIndels + 5 ; i++) {
            indels[i-5] = sample.getCounts(i).getGenotypeCountForwardStrand() + sample.getCounts(i).getGenotypeCountReverseStrand();
        }
        return String.format("mutated (%s) sample counts A=%d T=%d C=%d G=%d N=%d %s indels:%s",mutatedAllele,a,t,c,g,n,fb, Arrays.toString(indels));
    }


}