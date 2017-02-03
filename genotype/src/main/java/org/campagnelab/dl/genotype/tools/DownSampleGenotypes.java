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
 * Down-sample genotypes that do not have certain characteristics. Used to over-sample SNPs and het variants.
 *
 */
public class DownSampleGenotypes extends AbstractTool<DownSampleGenotypeArguments> {


    private RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    private Random random;
    private int recordsIncluded = 0;
    private int sitesNotSampled = 0;
    private int inputNumRecords = 0;

    static private Logger LOG = LoggerFactory.getLogger(DownSampleGenotypes.class);
    private int numIndel;
    private int numHeterozygotes;


    public static void main(String[] args) {

        DownSampleGenotypes tool = new DownSampleGenotypes();
        tool.parseArguments(args, "DownSampleGenotypes", tool.createArguments());
        tool.execute();
    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {
        // Check that the arguments are sound:
        if (args().keepAllHeterozygotes) {
            System.out.println("Will keep heterozygous sites.");
        }
        if (args().keepAllIndels) {
            System.out.println("Will keep indel sites.");
        }
        if (args().otherSamplingRate == 0 && !args().keepAllHeterozygotes && !args().keepAllIndels) {
            System.out.println("These arguments would result in nothing written to the output. Aborting.");
            System.exit(1);
        }

        //init random
        random = new XorShift1024StarRandom(args().seed);
        //get reference genome
        String genomePath = args().genomeFilename;
        try {
            System.err.println("Loading genome cache " + genomePath);
            genome = new RandomAccessSequenceCache();
            genome.load(genomePath, "min", "max");
            System.err.println("Done loading genome. ");
        } catch (Exception e) {
            throw new RuntimeException("Could not load genome cache " + args().genomeFilename, e);
        }

        try {
            int numSkippedSoFar = 0;
            int numKeptSoFar = 0;
            int numOtherKeptSoFar = 0;
            int numMustKeepSoFar = 0;
            RecordReader source = new RecordReader(args().inputFile);
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            recordLogger.start();

            for (BaseInformationRecords.BaseInformation rec : source) {
                int referenceIndex = genome.getReferenceIndex(rec.getReferenceId());
                char referenceBase = genome.get(referenceIndex, rec.getPosition());
                String refBase = Character.toString(referenceBase);
                inputNumRecords++;
                boolean isVariant = rec.getSamples(0).getIsVariant();
                String trueGenotype = rec.getTrueGenotype();

                boolean keep = false;

                //todo verify use of isIndel here
                boolean indel = GenotypeHelper.isIndel(refBase,trueGenotype);
                boolean heterozygote = GenotypeHelper.isHeterozygote(trueGenotype);

                keep |= args().keepAllIndels && indel;
                keep |= args().keepAllHeterozygotes && heterozygote;
                if (keep) {
                    numMustKeepSoFar++;
                }
                if (!keep) {

                    if (random.nextFloat() > args().otherSamplingRate) {
                        sitesNotSampled++;
                        numOtherKeptSoFar++;
                    } else {
                        keep = true;

                    }
                    if (args().balancingRatio != null) {
                        /**
                        System.out.printf("sitesNotSampled %d numKeptSoFar %d otherSamplingRate %f%n",sitesNotSampled,numKeptSoFar, args().otherSamplingRate);

                        float hasBeen=  numKeptSoFar-numMustKeepSoFar;
                        float shouldBe = numKeptSoFar*args().balancingRatio;
                        System.out.printf("hasBeen=%f shouldBe=%f%n",hasBeen,shouldBe);
                        args().otherSamplingRate = Math.max(0, hasBeen  / shouldBe );
                   */
                        throw new UnsupportedOperationException("This argument is currently not supported.");
                    }
                }
                if (keep) {
                    if (indel) {
                        numIndel++;
                    }
                    if (heterozygote) {
                        numHeterozygotes++;
                    }
                    recordsIncluded++;
                    dest.appendEntry(rec);
                    numKeptSoFar++;
                } else {
                    numSkippedSoFar++;
                }
                recordLogger.update();
            }
            recordLogger.done();
            dest.setCustomProperties(getStatProperties(source));
            dest.close();
            printStats(numSkippedSoFar);
        } catch (IOException e) {
            System.err.println("IO exception, perhaps sbi file not found?");
            e.printStackTrace();
            System.exit(1);
        }
    }


    @Override
    public DownSampleGenotypeArguments createArguments() {
        return new DownSampleGenotypeArguments();
    }


    private void printStats(int numSkippedSoFar) {
        System.out.println(numSkippedSoFar + " number of sites removed from the file..");
        System.out.println(recordsIncluded + " labeled records written.");
        System.out.println(numHeterozygotes + " heterozygotes records written.");
        System.out.println(numIndel + " indels records written.");
    }


    public Properties getStatProperties(RecordReader source) {
        Properties result = source.getProperties();
        result.put("downSampleGenotypes.sitesNotSampled", Integer.toString(sitesNotSampled));
        result.put("downSampleGenotypes.numIndelsWritten", Integer.toString(numIndel));
        result.put("downSampleGenotypes.numHeterezygotesWritten", Integer.toString(numHeterozygotes));
        result.put("downSampleGenotypes.sitesNotSampled", Integer.toString(sitesNotSampled));
        result.put("downSampleGenotypes.input.numRecords", Integer.toString(inputNumRecords));
        result.put("downSampleGenotypes.otherSamplingRate", Float.toString(args().otherSamplingRate));
        return result;
    }
}
