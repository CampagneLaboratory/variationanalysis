package org.campagnelab.dl.varanalysis.tools;

import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Random;

/**
 * Split a BSI file into several parts. Useful for creating training/validation/test splits of a large dataset.split
 * Created by fac2003 on 9/2/16.
 */
public class Split extends AbstractTool<SplitArguments>{


    static private Logger LOG = LoggerFactory.getLogger(Split.class);

    public static void main(String[] args) {

        Split tool = new Split();
        tool.parseArguments(args, "Split", tool.createArguments());
        tool.execute();
    }
    @Override
    public void execute()  {
        if (args().fractions.size() != args().suffixes.size()) {
            System.err.println("You must give exactly the same number of fractions and suffixes.");
            System.exit(1);
        }

        int numOutputs = args().fractions.size();
        if (numOutputs == 1) {
            System.err.println("Splitting a file into one fraction is not useful. Aborting.");
            System.exit(1);
        }
        RecordWriter outputWriters[] = new RecordWriter[numOutputs];

        try (RecordReader reader = new RecordReader(args().inputFile)) {
            double[] fractions = new double[numOutputs];
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

            Random rand = new XoRoShiRo128PlusRandom();
//set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "records";
            pgRead.expectedUpdates = Math.min(reader.getTotalRecords(), args().writeN);
            pgRead.displayFreeMemory = true;
            pgRead.start();
            long numWritten = 0;
            for (BaseInformationRecords.BaseInformation record : reader) {
                double choice = rand.nextDouble();
                double cumulativeFration = 0;
                int index;
                for (index = 0; index < numOutputs; index++) {
                    cumulativeFration += fractions[index];
                    if (choice < cumulativeFration) {
                        // we found the index of the output that should contain this record.
                        break;
                    }
                }
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

   

    @Override
    public SplitArguments createArguments() {
            return new SplitArguments();
    }


}
