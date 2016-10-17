package org.campagnelab.dl.varanalysis.intermediaries;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * The randomizer object iterates over a parquet file and randomizes the order of records in batches.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class Randomizer2 extends Intermediary {

    static private Logger LOG = LoggerFactory.getLogger(Randomizer2.class);
    private List<RecordReader> sources = new ObjectArrayList<RecordReader>();
    private Randomizer2Arguments arguments;

    public static void main(String[] args) throws IOException {
        Randomizer2Arguments arguments = parseArguments(args, "Randomizer2");

        //randomize
        Randomizer2 r = new Randomizer2();
        r.arguments = arguments;
        r.setSourceFiles(arguments.inputFiles.toArray(new String[0]));
        r.executeOver(null, arguments.outputFile);

    }

    protected static Randomizer2Arguments parseArguments(String[] args, String commandName) {
        Randomizer2Arguments arguments = new Randomizer2Arguments();
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

    public void execute(String inPath, String outPath, int blockSize, int pageSize) throws IOException {
        String workingDir = new File(outPath).getParent();
        if (workingDir == null) {
            workingDir = ".";
        }
        try {
            int totalRecords = 0;
            for (String filename : sourceFilenames) {
                RecordReader source = new RecordReader(filename);
                totalRecords += source.getTotalRecords();
                source.close();
            }
            int numBuckets = (totalRecords / arguments.recordsPerBucket) + 1;

            new File(workingDir + "/tmp").mkdir();
            List<RecordWriter> bucketWriters = new ObjectArrayList<RecordWriter>(numBuckets);
            for (int i = 0; i < numBuckets; i++) {
                bucketWriters.add(new RecordWriter(workingDir + "/tmp/bucket" + i, arguments.chunkSizePerWriter));
            }
            RecordWriter allWriter = new RecordWriter(outPath);
            Random rand = new XoRoShiRo128PlusRandom();

            //set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "read";
            pgRead.expectedUpdates = totalRecords;
            pgRead.displayFreeMemory = true;
            pgRead.start();

            //fill buckets randomly
            System.out.println("Filling " + numBuckets + " temp buckets randomly");
            for (String filename : sourceFilenames) {
                RecordReader source = new RecordReader(filename);
                for (BaseInformationRecords.BaseInformation rec : source) {
                    int bucket = rand.nextInt(numBuckets);
                    bucketWriters.get(bucket).writeRecord(rec);
                    pgRead.lightUpdate();
                }
                source.close();
                System.gc();
            }

            pgRead.stop();

            System.out.println("Shuffling contents of each bucket and writing to output file");
            //iterate over buckets
            ProgressLogger pgTempBucket = new ProgressLogger(LOG);
            pgTempBucket.itemsName = "read";
            pgTempBucket.expectedUpdates = numBuckets;
            pgTempBucket.displayFreeMemory = true;
            pgTempBucket.start();
            int i = 0;
            for (RecordWriter bucketWriter : bucketWriters) {
                bucketWriter.close();

                //put contents of bucket in a list
                RecordReader bucketReader = new RecordReader(workingDir + "/tmp/bucket" + i );
                List<BaseInformationRecords.BaseInformation> records = new ObjectArrayList<>(arguments.recordsPerBucket);
                for (BaseInformationRecords.BaseInformation rec : bucketReader) {
                    records.add(rec);
                }
                bucketReader.close();

                //shuffle list
                Collections.shuffle(records, rand);

                //write list to final file
                for (BaseInformationRecords.BaseInformation rec : records) {
                    allWriter.writeRecord(rec);
                }
                i++;
                pgTempBucket.update();
            }
            pgTempBucket.stop();
            allWriter.close();

            //delete temp files
            FileUtils.deleteDirectory(new File((workingDir + "/tmp")));

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    String[] sourceFilenames;

    void setSourceFiles(String[] args) throws IOException {
        this.sourceFilenames = args;
    }

}
