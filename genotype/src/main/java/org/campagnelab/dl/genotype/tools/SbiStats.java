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
import java.util.Set;

/**
 * Use sbistats to see the distibution of types of records in a dataset.
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
                try {
            RecordReader source = new RecordReader(args().inputFile);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            int recordsLabeled = 0;
            recordLogger.start();

            int numIndels = 0;
            int numSnps = 0;
            int numHetSnps = 0;
            int numHomSnps = 0;
            int numHetIndels = 0;
            int numHomIndels = 0;
            int numVariants = 0;
            int numHasIndel = 0;

            for (BaseInformationRecords.BaseInformation rec : source) {
                boolean isIndel = false;
                boolean hasIndel = false;
                boolean isVariant = false;
                boolean isSnp = false;
                boolean heterozygous = false;
                for (BaseInformationRecords.CountInfo c : rec.getSamples(args().sampleIndex).getCountsList()){
                    if (c.getIsIndel() && c.getIsCalled()) {
                        isIndel = true;
                    }
                    if (c.getIsIndel()){
                        hasIndel = true;
                    }
                    if (!c.getToSequence().equals(c.getFromSequence())){
                        isVariant = true;
                    }
                    if (isVariant && !isIndel){
                        isSnp = true;
                    }
                }
                String trueGenotype = rec.getTrueGenotype();
                Set<String> alleles = GenotypeHelper.getAlleles(trueGenotype);
                if (alleles.size() > 1) {
                    heterozygous = true;
                }
                if (isVariant){
                    numVariants++;
                }
                if (isIndel){
                    numIndels++;
                    if(heterozygous){
                        numHetIndels++;
                    } else {
                        numHomIndels++;
                    }
                } if (isSnp){
                    numSnps++;
                    if(heterozygous){
                        numHetSnps++;
                    } else {
                        numHomSnps++;
                    }
                }
                if (hasIndel){
                    numHasIndel++;
                }
                recordLogger.lightUpdate();
            }
            recordLogger.done();
            System.out.println("numIndels = " + numIndels );
            System.out.println("numSnps = " + numSnps );
            System.out.println("numHetSnps = " + numHetSnps );
            System.out.println("numHomSnps = " + numHomSnps );
            System.out.println("numHetIndels = " + numHetIndels );
            System.out.println("numHomIndels = " + numHomIndels );
            System.out.println("numVariants = " + numVariants );
                    System.out.println("numHasIndel = " + numHasIndel );
            } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public SbiStatsArguments createArguments() {
        return new SbiStatsArguments();
    }


}
