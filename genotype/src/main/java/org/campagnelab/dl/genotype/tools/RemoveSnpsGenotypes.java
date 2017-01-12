package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.logging.ProgressLogger;
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
import java.util.Properties;
import java.util.Random;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class RemoveSnpsGenotypes extends AbstractTool<RemoveSnpsGenotypesArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    Random random;
    int recordsIncluded = 0;
    int snpsRemoved = 0;
    int referenceNotSampled = 0;
    int inputNumRecords = 0;
    float referenceSamplingRate;

    static private Logger LOG = LoggerFactory.getLogger(RemoveSnpsGenotypes.class);

    public static void main(String[] args) {

        RemoveSnpsGenotypes tool = new RemoveSnpsGenotypes();
        tool.parseArguments(args, "RemoveSnpsGenotypes", tool.createArguments());
        tool.execute();


    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {

        //init random
        random = new Random(args().seed);
        referenceSamplingRate = args().referenceSamplingRate;

        try {
            System.err.println("Done loading genome. ");
            RecordReader source = new RecordReader(args().inputFile);
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            recordLogger.start();
            for (BaseInformationRecords.BaseInformation rec : source) {
                inputNumRecords++;
                boolean isVariant = rec.getSamples(0).getIsVariant();
                boolean isSnp = isVariant && !GenotypeHelper.isIndel(rec.getTrueGenotype());
                if (isSnp){
                    snpsRemoved++;
                    continue;
                }
                if (!isVariant) {
                    if (random.nextFloat() > referenceSamplingRate ) {
                        referenceNotSampled++;
                        continue;
                    }
                }
                recordsIncluded++;
                dest.appendEntry(rec);
                recordLogger.lightUpdate();
            }
            recordLogger.done();
            dest.setCustomProperties(getStatProperties(source));
            dest.close();
            printStats(recordsIncluded,snpsRemoved,referenceNotSampled);
        } catch (IOException e) {
            System.err.println("IO exception, perhaps sbi file not found?");
            e.printStackTrace();
            System.exit(1);
        }
    }


    @Override
    public RemoveSnpsGenotypesArguments createArguments() {
        return new RemoveSnpsGenotypesArguments();
    }


    private void printStats(int recordsIncluded, int snpsRemoved, int referenceNotSampled) {
        System.out.println(referenceNotSampled + " number of reference examples not sampled/written to file");
        System.out.println(snpsRemoved + " number of snps removed from the file..");
        System.out.println(recordsIncluded + " labeled records written.");

    }


    public Properties getStatProperties(RecordReader source) {
        Properties result = source.getProperties();
        result.put("removeSnps.numSnpsRemoved", Integer.toString(snpsRemoved));
        result.put("removeSnps.referenceNotSampled", Integer.toString(referenceNotSampled));
        result.put("removeSnps.input.numRecords", inputNumRecords);
        result.put("removeSnps.referenceSamplingRate", Float.toString(referenceSamplingRate));
        return result;
    }
}
