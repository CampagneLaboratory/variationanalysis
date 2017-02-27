package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.fastutil.doubles.DoubleArraySet;
import it.unimi.dsi.fastutil.doubles.DoubleSet;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Random;

/**
 * Split a BSI file into several parts. Useful for creating training/validation/test splits of a large dataset.split
 * Created by fac2003 on 9/2/16.
 */
public class Split extends AbstractTool<SplitArguments> {


    static private Logger LOG = LoggerFactory.getLogger(Split.class);
    private int numOutputs;
    private double[] fractions;
    private DoubleSet excludedIndices = new DoubleArraySet();

    public static void main(String[] args) {

        Split tool = new Split();
        tool.parseArguments(args, "Split", tool.createArguments());
        tool.execute();
    }

    Random rand;
    Object2IntMap<String> chomosomeToSuffixIndex;

    @Override
    public void execute() {
        if (args().fractions.size() != args().suffixes.size()) {
            System.err.println("You must give exactly the same number of fractions and suffixes.");
            System.exit(1);
        }
        chomosomeToSuffixIndex = new Object2IntArrayMap<>();
        chomosomeToSuffixIndex.defaultReturnValue(-1);
        if (args().destinationOverride != null) {
            // parsing destinationOverride:
            String[] tokens = args().destinationOverride.split(":");
            if (tokens.length != 2) {
                System.err.println("Invalid syntax for destinationOverride. Must have two tokens separated by colon: suffix and list of chromosomes, such as test:chr1,chr2");
                System.exit(1);
            }
            String suffix = tokens[0];
            int indexOfSuffix = args().suffixes.indexOf(suffix);
            if (indexOfSuffix == -1) {
                System.err.printf("Invalid syntax for destinationOverride: Suffix %s is unknown.", suffix);
                System.exit(1);
            }
            String[] chromosomes = tokens[1].split(",");
            for (String chromosome : chromosomes) {
                chomosomeToSuffixIndex.put(chromosome, indexOfSuffix);
                // exclude the chromosome from regular sampling:
                excludedIndices.add(indexOfSuffix);
            }
        }
        numOutputs = args().fractions.size();
        if (numOutputs == 1) {
            System.err.println("Splitting a file into one fraction is not useful. Aborting.");
            System.exit(1);
        }
        RecordWriter outputWriters[] = new RecordWriter[numOutputs];

        try (RecordReader reader = new RecordReader(args().inputFile)) {
            fractions = new double[numOutputs];
            double sumFractions = 0;
            for (int i = 0; i < numOutputs; i++) {
                outputWriters[i] = new RecordWriter(args().outputFile + args().suffixes.get(i));
                fractions[i] = args().fractions.get(i);
                sumFractions += fractions[i];
            }
            // normalize fractions:
            for (int i = 0; i < numOutputs; i++) {
                fractions[i] /= sumFractions;
            }

            rand = new XoRoShiRo128PlusRandom();
//set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "records";
            pgRead.expectedUpdates = Math.min(reader.getTotalRecords(), args().writeN);
            pgRead.displayFreeMemory = true;
            pgRead.start();
            long numWritten = 0;
            for (BaseInformationRecords.BaseInformation record : reader) {
                int index = recorgBelongsTo(record);

                outputWriters[index].writeRecord(record);
                pgRead.lightUpdate();
                numWritten += 1;
                if (numWritten > args().writeN) {
                    break;
                }
            }
            pgRead.stop();
            for (int i = 0; i < numOutputs; i++) {
                outputWriters[i].close();
            }
            reader.close();
        } catch (IOException e) {
            System.err.println("Unable to load or write files. Check command line arguments.");
        }
    }


    private int recorgBelongsTo(BaseInformationRecords.BaseInformation record) {
        final String chromosome = record.getReferenceId();
        //      System.out.println(chromosome);

        int overrideIndex = chromosome == null ? -1 : chomosomeToSuffixIndex.getInt(chromosome);
        if (overrideIndex == -1 || args().destinationOverride == null) {

            double choice = rand.nextDouble();
            double cumulativeFration = 0;
            int index;
            for (index = 0; index < numOutputs; index++) {
                cumulativeFration += fractions[index];
                if (choice < cumulativeFration) {
                    if (excludedIndices.contains(index)) {
                        index = 0;
                        continue;
                    }
                    // we found the index of the output that should contain this record.
                    break;
                }
            }
            return index;
        } else {
            // System.out.println("Overriding destination: " + overrideIndex);
            return overrideIndex;
        }
    }


    @Override
    public SplitArguments createArguments() {
        return new SplitArguments();
    }


}
