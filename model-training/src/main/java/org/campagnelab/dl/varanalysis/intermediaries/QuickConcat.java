package org.campagnelab.dl.varanalysis.intermediaries;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.compression.MessageChunksWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;

/**
 * A utility to quickly concatenate a list of .sbi files. Concatenation does not decompress each file and simply
 * concatenates the bytes and add up the number of records.
 */
public class QuickConcat {
    private QuickConcatArguments arguments;
    static private Logger LOG = LoggerFactory.getLogger(QuickConcat.class);

    public static void main(String[] args) throws IOException {
        QuickConcatArguments arguments = parseArguments(args, "QuickConcat");
        QuickConcat r = new QuickConcat();
        r.arguments = arguments;

        r.performQuickConcat(arguments.inputFiles.toArray(new String[0]),arguments.outputFile);

    }

    protected static QuickConcatArguments parseArguments(String[] args, String commandName) {
        QuickConcatArguments arguments = new QuickConcatArguments();
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
    /**
     * This version does a quick concat. It does NO filtering. It gathers no stats,
     * but, will quickly concat multiple compact-reads files together using NIO.
     * It should be noted that this method is >MUCH< faster.
     * Copy all of the input files except the last MessageChunksWriter.DELIMITER_LENGTH
     * bytes of the first n-1 input files and the entire last input file
     * to the output file.
     * @throws IOException
     * @param inputFilenames
     * @param outputBasename
     */
    private void performQuickConcat(String[] inputFilenames, String outputBasename) throws IOException {
        System.out.println("quick concatenating files");
        File outputFile = new File(outputBasename);
        if (outputFile.exists()) {
            System.err.println("The output file already exists. Please delete it before running concat.");
            return;
        }
        outputFile.createNewFile();
        //set up logger
        ProgressLogger progressLogger = new ProgressLogger(LOG);
        progressLogger.itemsName = "files";
        progressLogger.expectedUpdates = inputFilenames.length;
        progressLogger.displayFreeMemory = true;
        progressLogger.start();
        FileChannel input = null;
        FileChannel output = null;
        long bufferSize = arguments.copyBufferSize;

        long totalElements=0;
        for (final String inputFilename : inputFilenames) {
            SequenceBaseInformationReader reader=new SequenceBaseInformationReader(inputFilename);
            totalElements+=  reader.getTotalRecords();
            reader.close();
        }
        SequenceBaseInformationWriter.writeProperties(outputBasename, totalElements);
        try {
            output = new FileOutputStream(outputFile+".sbi").getChannel();
            int lastFileNumToCopy = inputFilenames.length - 1;
            int curFileNum = 0;
            for (final String inputFilename : inputFilenames) {
                System.out.printf("Reading from %s%n", inputFilename);
                input = new FileInputStream(inputFilename).getChannel();
                long bytesToCopy = input.size();
                if (curFileNum++ < lastFileNumToCopy) {
                    // Compact-reads files end with a delimiter (8 x 0xff)
                    // followed by a 4 byte int 0 (4 x 0x00). Strip
                    // these on all but the last file.
                    bytesToCopy -= (MessageChunksWriter.DELIMITER_LENGTH + 1+ MessageChunksWriter.SIZE_OF_MESSAGE_LENGTH);
                }

                // Copy the file about 10 megabytes at a time. It would probably
                // be marginally faster to just tell NIO to copy the ENTIRE file
                // in one go, but with very large files Java will freeze until the
                // entire chunck is copied so this makes for a more responsive program
                // should you want to ^C in the middle of the copy. Also, with the single
                // transferTo() you might not see any file size changes in the output file
                // until the entire copy is complete.
                long position = 0;
                while (position < bytesToCopy) {
                    long bytesToCopyThisTime = Math.min(bufferSize, bytesToCopy - position);
                    position += input.transferTo(position, bytesToCopyThisTime, output);
                }
                input.close();
                input = null;
                progressLogger.update();
            }
            System.out.printf("Concatenated %d files.%n", lastFileNumToCopy + 1);
            progressLogger.stop();
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
            }
        }
    }
}
