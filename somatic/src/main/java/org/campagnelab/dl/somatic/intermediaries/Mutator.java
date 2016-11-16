package org.campagnelab.dl.somatic.intermediaries;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * The mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class Mutator extends Intermediary {
    static private Logger LOG = LoggerFactory.getLogger(Mutator.class);


    final boolean MUTATE = true;


    public static void main(String[] args) throws IOException {
        //new ParquetPrinter(args[0]).print();
        if (args.length < 3) {
            System.err.println("usage: input.parquet mutated-filename mutated-randomized-filename");
            System.exit(1);
        }
        new Randomizer().executeOver(args[0], args[1]);
        new Mutator().executeOver(args[1], args[2]);
        //new ParquetPrinter(args[2]).print();
    }

    public Mutator() {
        strategy = new FirstSimulationStrategy(2323);
    }

    public void setSeed(int seed) {

        strategy.setSeed(seed);
    }

    public Mutator(int deltaSmall, int deltaBig, int zygHeuristic, long seed) {
        strategy = new FirstSimulationStrategy(deltaSmall, deltaBig, zygHeuristic,seed);
        strategy.setSeed(2323);
    }

    public void execute(String in, String out, int blockSize, int pageSize) throws IOException {
        RecordReader reader = new RecordReader(in);
        RecordWriter writer = new RecordWriter(out);

        //set up logger
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = "mutation";
        pgReadWrite.expectedUpdates = reader.getTotalRecords();
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        for (BaseInformationRecords.BaseInformation base : reader) {
            BaseInformationRecords.BaseInformation.Builder baseBuilder = base.toBuilder();

            //set is tumor flags (false sample 0 (germline), true sample 1 (somatic)
            baseBuilder.setSamples(0, baseBuilder.getSamples(0).toBuilder().setIsTumor(false));
            baseBuilder.setSamples(1, baseBuilder.getSamples(1).toBuilder().setIsTumor(true));
            writer.writeRecord(baseBuilder.build());

            //mutate record and write it again
            if (MUTATE) {
                writer.writeRecord(mutate(baseBuilder));
            }
            pgReadWrite.update();
        }
        pgReadWrite.stop();
        reader.close();
        writer.close();

    }

    SimulationStrategy strategy ;

    //backward strand appended to forward strand in input and output
    //10 fields: ATCGOATCGO
    //List instead of array because of avro code generation...
    protected BaseInformationRecords.BaseInformation mutate(BaseInformationRecords.BaseInformation record) {

        return strategy.mutate(true, record, record.getSamples(0), record.getSamples(1), null);
    }


    public BaseInformationRecords.BaseInformation mutate(BaseInformationRecords.BaseInformation.Builder builder) {
        return mutate(builder.build());
    }
}