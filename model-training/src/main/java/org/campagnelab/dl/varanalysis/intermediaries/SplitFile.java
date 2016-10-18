package org.campagnelab.dl.varanalysis.intermediaries;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
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
public class SplitFile {


    static private Logger LOG = LoggerFactory.getLogger(SplitFile.class);

    public static void main(String[] args) throws IOException {
        SplitFileArguments arguments = parseArguments(args, "SplitFile");

        SplitFile splitter = new SplitFile();
        splitter.execute(arguments);
    }

    private void execute(SplitFileArguments arguments) throws IOException {
        if (arguments.fractions.size() != arguments.suffixes.size()) {
            System.err.println("You must give exactly the same number of fractions and suffixes.");
            System.exit(1);
        }

        int numOutputs = arguments.fractions.size();
        if (numOutputs == 1) {
            System.err.println("Splitting a file into one fraction is not useful. Aborting.");
            System.exit(1);
        }
        RecordWriter outputWriters[] = new RecordWriter[numOutputs];

        RecordReader reader = new RecordReader(arguments.inputFile);
        double[] fractions = new double[numOutputs];
        double sumFractions = 0;
        for (int i = 0; i < numOutputs; i++) {
            outputWriters[i] = new RecordWriter(arguments.outputFile + arguments.suffixes.get(i));
            fractions[i] = arguments.fractions.get(i);
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
        pgRead.expectedUpdates = Math.min(reader.getTotalRecords(),arguments.writeN);
        pgRead.displayFreeMemory = true;
        pgRead.start();
        long numWritten=0;
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
            numWritten+=1;
            if (numWritten>arguments.writeN) {
                break;
            }
        }
        pgRead.stop();
        for (int i = 0; i < numOutputs; i++) {
            outputWriters[i].close();
        }
        reader.close();
    }

    protected static SplitFileArguments parseArguments(String[] args, String commandName) {
        SplitFileArguments arguments = new SplitFileArguments();
        JCommander commander = new JCommander(arguments);
        commander.setProgramName(commandName);
        try {
            commander.parse(args);
        } catch (ParameterException e) {

            commander.usage();
            throw e;

        }
        return arguments;
    }
}
