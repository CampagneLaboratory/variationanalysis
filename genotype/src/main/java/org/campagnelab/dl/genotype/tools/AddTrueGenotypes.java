package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.AddTrueGenotypeHelper;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class AddTrueGenotypes extends AbstractTool<AddTrueGenotypesArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    final public static boolean PRINT_INDEL_ERROR_CONTEXT = false;

    static private Logger LOG = LoggerFactory.getLogger(AddTrueGenotypes.class);

    public static void main(String[] args) {

        AddTrueGenotypes tool = new AddTrueGenotypes();
        tool.parseArguments(args, "AddTrueGenotypes", tool.createArguments());
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
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);
            AddTrueGenotypeHelper addTrueGenotypeHelper = new AddTrueGenotypeHelper();
            addTrueGenotypeHelper.configure(
                    args().genotypeMap,
                    genome,
                    args().sampleIndex,
                    args().considerIndels,
                    args().indelsAsRef,
                    args().referenceSamplingRate);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            int recordsLabeled = 0;
            recordLogger.start();
            ObjectArrayList<BaseInformationRecords.BaseInformation> recContext = new ObjectArrayList<>(1000);
            for (BaseInformationRecords.BaseInformation rec : source) {
                boolean keep = false;
                if (PRINT_INDEL_ERROR_CONTEXT){
                    recContext.add(rec);
                    if (recContext.size() < 50){
                        continue;
                    }
                    keep = addTrueGenotypeHelper.addTrueGenotype(recContext.get(recContext.size()/2),recContext);
                    if (keep) {
                        dest.appendEntry(addTrueGenotypeHelper.labeledEntry());
                    }

                    recContext.remove(0);
                } else {
                    keep = addTrueGenotypeHelper.addTrueGenotype(rec);
                    if (keep) {
                        dest.appendEntry(addTrueGenotypeHelper.labeledEntry());
                    }
                }
                recordsLabeled++;
                recordLogger.lightUpdate();


            }
            recordLogger.done();
            dest.setCustomProperties(addTrueGenotypeHelper.getStatProperties());
            dest.close();
            addTrueGenotypeHelper.printStats();
            } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public AddTrueGenotypesArguments createArguments() {
        return new AddTrueGenotypesArguments();
    }


}
