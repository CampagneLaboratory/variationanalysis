package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
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
 * Filter an SBI dataset to keep sites that match certain criteria.
 * <p>
 *
 *
 * @author Fabien Campagne
 */
public class FilterSBI extends AbstractTool<FilterSBIArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    Random random;
    long filteredCount=0;
    long includedCount=0;
    double otherSamplingRate;


    static private Logger LOG = LoggerFactory.getLogger(FilterSBI.class);

    public static void main(String[] args) {

        FilterSBI tool = new FilterSBI();
        tool.parseArguments(args, "FilterSBI", tool.createArguments());
        tool.execute();


    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {

        //init random
        random = new XorShift1024StarRandom(args().seed);
        double otherSamplingRate = args().otherSamplingRate;

        try {

            RecordReader source = new RecordReader(args().inputFile);
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to filter");
            recordLogger.start();
            for (BaseInformationRecords.BaseInformation rec : source) {

                boolean isVariant = rec.getSamples(0).getIsVariant();
                boolean isSnp = isVariant && !GenotypeHelper.isIndel(rec.getReferenceBase(), rec.getTrueGenotype());
                boolean toRemove = false;
                if (args().removeSNPs && isSnp) {
                    toRemove = true;
                    filteredCount++;
                }
                if (args().removeReferenceMatching && !isVariant) {

                    toRemove = true;
                }

                if (toRemove) {
                    if (random.nextFloat() < otherSamplingRate) {
                        //override the remove flag. We want to keep a subset of those we remove.
                        toRemove = false;
                        filteredCount--;
                    }
                }
                if (!toRemove) {
                    includedCount++;
                    dest.appendEntry(rec);
                }
                recordLogger.lightUpdate();
            }
            recordLogger.done();
            dest.setCustomProperties(getStatProperties(source));
            dest.close();
            printStats(includedCount, filteredCount);
        } catch (IOException e) {
            System.err.println("IO exception, perhaps the sbi file was not found?");
            e.printStackTrace();
            System.exit(1);
        }
    }


    @Override
    public FilterSBIArguments createArguments() {
        return new FilterSBIArguments();
    }


    private void printStats(long recordsIncluded, long filteredCount) {
        System.out.printf( "included: %d%n",recordsIncluded);
        System.out.printf(  "sites filtered: %d%n",filteredCount);
        System.out.printf(  "Removed %f %% of the input sites.%n",100*(0.0+filteredCount)/(0.0+filteredCount+includedCount));

    }


    public Properties getStatProperties(RecordReader source) {
        Properties result = source.getProperties();
        result.put("removeSnps.filtered", Long.toString(filteredCount));
        result.put("removeSnps.included", Long.toString(includedCount));

        return result;
    }
}
