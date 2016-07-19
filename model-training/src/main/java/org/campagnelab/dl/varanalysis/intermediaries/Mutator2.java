package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift128PlusRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

/**
 * The mutator object iterates over a parquet file and creates additional copies of every record, where
 * some copies are mutated, and some are not. In constrast to Mutator, Mutator2 observes the sum of counts in the
 * second sample, discards the counts in the second sample, and recreates them using a simulation protocol
 * that either simulates random sampling from the same genotype, or from a second tumor genotype derived from germline.
 * <p>
 *
 * @author Fabien Campagne
 */
public class Mutator2 extends Intermediary {
    private static final int CHUNK_SIZE = 10000;
    private static final int NUM_SIMULATED_RECORD_PER_DATUM = 10;
    static private Logger LOG = LoggerFactory.getLogger(Mutator2.class);

    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record (ceilinged fraction)
    double deltaSmall = 0.0;
    //max fraction of bases mutated at a record (floored fraction)
    double deltaBig = 1;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;
    final String[] STRING = new String[]{"A", "T", "C", "G"};
    Random rand;
    final boolean MUTATE = true;
    SimulationStrategy stragegy = new SimulationStrategyImpl();

    public static void main(String[] args) throws IOException {
        //new ParquetPrinter(args[0]).print();
        if (args.length < 3) {
            System.err.println("usage: input.parquet mutated-filename mutated-randomized-filename");
            System.exit(1);
        }
        new Randomizer().executeOver(args[0], args[1]);
        new Mutator2().executeOver(args[1], args[2]);
        //new ParquetPrinter(args[2]).print();
    }

    public Mutator2() {
        setSeed(2323);
    }

    public void setSeed(int seed) {

        rand = new XorShift128PlusRandom(seed);
    }

    public Mutator2(int deltaSmall, int deltaBig, int zygHeuristic) {
        this.deltaSmall = deltaSmall;
        this.deltaBig = deltaBig;
        this.zygHeuristic = zygHeuristic;
        setSeed(2323);
    }

    public void execute(String in, String out, int blockSize, int pageSize) throws IOException {
        RecordReader reader = new RecordReader(in);
        RecordWriter writer = new RecordWriter(out, blockSize, pageSize, true);

        //set up logger
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = "mutation";
        pgReadWrite.expectedUpdates = reader.getTotalRecords();
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();
        SimulationCharacteristics sim = new SimulationCharacteristics();
        try {
            for (BaseInformationRecords.BaseInformation base : reader) {

                sim.observe(base);
                if (sim.size() >= CHUNK_SIZE) {
                    sim.batchIsComplete();
                    processBatch(sim, writer);
                    sim.clear();
                }

                pgReadWrite.lightUpdate();
            }
        } finally {
            pgReadWrite.stop();
            reader.close();
            writer.close();

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

            for (int i = 0; i < NUM_SIMULATED_RECORD_PER_DATUM; i++) {

                shufflingList.add(stragegy.mutate(true, record, record.getSamples(0), record.getSamples(1), sim));
                shufflingList.add(stragegy.mutate(false, record, record.getSamples(1), record.getSamples(0), sim));
            }
        }
        Collections.shuffle(shufflingList);
        for (BaseInformationRecords.BaseInformation record : shufflingList) {
            writer.writeRecord(record);
        }

    }


}