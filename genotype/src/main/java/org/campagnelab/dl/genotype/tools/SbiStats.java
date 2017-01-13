package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.AddTrueGenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class SbiStats extends AbstractTool<SbiStatsArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();

    static private Logger LOG = LoggerFactory.getLogger(SbiStats.class);

    public static void main(String[] args) {

        SbiStats tool = new SbiStats();
        tool.parseArguments(args, "SbiStatsArguments", tool.createArguments());
        tool.execute();

    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {

        //get reference genome
        String genomePath = args().genomeFilename;
        try {
            System.err.println("Loading genome cache " + genomePath);
            genome = new RandomAccessSequenceCache();
            genome.load(genomePath, "min", "max");
            System.err.println("Done loading genome. ");
        } catch (ClassNotFoundException e) {
            System.err.println("Could not load genome cache");
            e.printStackTrace();
            System.exit(1);
        } catch (IOException e) {
            System.err.println("Could not load genome cache");
            e.printStackTrace();
            System.exit(1);
        }


        try {
            RecordReader source = new RecordReader(args().inputFile);
            AddTrueGenotypeHelper addTrueGenotypeHelper = new AddTrueGenotypeHelper();
            addTrueGenotypeHelper.configure(
                    args().genotypeMap,
                    genome,
                    args().sampleIndex,
                    args().considerIndels,
                    1.0f);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            int recordsLabeled = 0;
            recordLogger.start();
            for (BaseInformationRecords.BaseInformation rec : source) {
                addTrueGenotypeHelper.addTrueGenotype(rec);

                recordLogger.lightUpdate();
            }
            recordLogger.done();
            addTrueGenotypeHelper.printStats();
            } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public SbiStatsArguments createArguments() {
        return new SbiStatsArguments();
    }


}
