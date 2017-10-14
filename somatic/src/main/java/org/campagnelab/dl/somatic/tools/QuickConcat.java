package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.compress.utils.IOUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.intermediaries.QuickConcatArguments;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.compression.MessageChunksWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Properties;

/**
 * A utility to quickly concatenate a list of .sbi files. Concatenation does not decompress each file and simply
 * concatenates the bytes and add up the number of records.
 */
public class QuickConcat extends AbstractTool<QuickConcatArguments> {
    static private Logger LOG = LoggerFactory.getLogger(QuickConcat.class);

    @Override
    public QuickConcatArguments createArguments() {
        return new QuickConcatArguments();
    }

    @Override
    public void execute() {

        performQuickConcat(arguments.inputFiles.toArray(new String[0]), arguments.outputFile);

    }

    public static void main(String[] args) throws IOException {
        QuickConcat r = new QuickConcat();
        r.parseArguments(args, "QuickConcat", r.createArguments());
        r.execute();

    }


    /**
     * This version does a quick concat. It does NO filtering. It gathers no stats,
     * but, will quickly concat multiple files together using NIO.
     * It should be noted that this method is >MUCH< faster.
     * Copy all of the input files except the last MessageChunksWriter.DELIMITER_LENGTH
     * bytes of the first n-1 input files and the entire last input file
     * to the output file.
     *
     * @param inputFilenames
     * @param outputBasename
     * @throws IOException
     */
    private void performQuickConcat(String[] inputFilenames, String outputBasename) {
        if (outputBasename.endsWith(".sbi")) {
            outputBasename = FilenameUtils.removeExtension(outputBasename);
        }
        System.out.println("quick concatenating files");
        inputFilenames = reorderFilenames(inputFilenames);

        File outputFile = new File(outputBasename);
        if (!args().force) {
            if (outputFile.exists()) {
                System.err.println("The output file already exists. Please delete it before running concat.");
                return;
            }
            try {
                outputFile.createNewFile();
            } catch (IOException e) {
                throw new RuntimeException("Unable to create destination file", e);
            }
        }

        //set up logger
        ProgressLogger progressLogger = new ProgressLogger(LOG);
        progressLogger.itemsName = "files";
        progressLogger.expectedUpdates = inputFilenames.length;
        progressLogger.displayFreeMemory = true;
        progressLogger.start();
        FileChannel input = null;
        FileChannel output = null;
        long bufferSize = arguments.copyBufferSize;

        ObjectList<Properties> properties = new ObjectArrayList<>();
        for (final String inputFilename : inputFilenames) {
            try {
                SequenceBaseInformationReader reader = new SequenceBaseInformationReader(inputFilename);
                properties.add(reader.getProperties());
                reader.close();
            } catch (IOException e) {
                throw new RuntimeException("Unable to open " + inputFilename, e);
            }
        }

        try {
            SequenceBaseInformationWriter.writeProperties(outputBasename, properties);
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to write properties", e);
        }

        try {
            output = new FileOutputStream(outputFile + ".sbi").getChannel();
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
                    bytesToCopy -= (MessageChunksWriter.DELIMITER_LENGTH + 1 + MessageChunksWriter.SIZE_OF_MESSAGE_LENGTH);
                }

                // Copy the file about 10 megabytes at a time. It would probably
                // be marginally faster to jugit st tell NIO to copy the ENTIRE file
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
        } catch (Exception e) {
            throw new RuntimeException("Unable to concatenate", e);
        } finally {
            IOUtils.closeQuietly(input);
            IOUtils.closeQuietly(output);
        }
    }

    private String[] reorderFilenames(String[] inputFilenames) {
        Integer[] order = new Integer[inputFilenames.length];
        int index = 0;
        boolean hasOrder = false;
        for (String inputFilename : inputFilenames) {
            String[] tokens = inputFilename.split("[\\-]");
            if ((tokens.length >= 2)) {
                try {
                    String lastToken = tokens[tokens.length - 1];
                    lastToken = lastToken.replace(".sbi", "");
                    int value = Integer.parseInt(lastToken);
                    order[index] = value;
                    hasOrder |= true;
                } catch (NumberFormatException e) {
                    order[index] = null;
                    hasOrder = false;
                }
                index++;
            } else {
                hasOrder = false;
            }
        }
        if (hasOrder) {
            System.out.println("Filenames follow pattern *-int.sbi. Concat will be in the increasing order of the integer values.");
            Comparator<? super String> compareFilenames = new Comparator<String>() {
                @Override
                public int compare(String filename1, String filename2) {
                    String[] tokens1 = filename1.split("[\\-]");
                    String[] tokens2 = filename2.split("[\\-]");
                    String lastToken1 = tokens1[tokens1.length - 1];
                    lastToken1 = lastToken1.replace(".sbi", "");

                    String lastToken2 = tokens2[tokens2.length - 1];
                    lastToken2 = lastToken2.replace(".sbi", "");

                    int value1 = Integer.parseInt(lastToken1);
                    int value2 = Integer.parseInt(lastToken2);
                    return Integer.compare(value1, value2);

                }
            };
            Arrays.sort(inputFilenames, compareFilenames);
            System.out.println("Sorted filename order:" + ObjectArrayList.wrap(inputFilenames).toString());
        } else {
            System.out.println("Filenames do not follow the pattern *-int.sbi, concat in the order provided on the command line.");
        }
        return inputFilenames;
    }


}
