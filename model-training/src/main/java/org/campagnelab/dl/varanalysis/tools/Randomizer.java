package org.campagnelab.dl.varanalysis.tools;


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
public class Randomizer extends AbstractTool<RandomizerArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Randomizer.class);

    public static void main(String[] args) {

        Randomizer tool = new Randomizer();
        tool.parseArguments(args, "Randomizer2", tool.createArguments());
        tool.execute();
    }


    @Override
    public void execute() {
        String workingDir = new File(args().outputFile).getParent();
        if (workingDir == null) {
            workingDir = ".";
        }
        try {
            long totalRecords = 0;
            for (String filename : args().inputFiles) {
                RecordReader source = new RecordReader(filename);
                totalRecords += source.getTotalRecords();
                source.close();
            }
            int numBuckets = (int)(totalRecords / arguments.recordsPerBucket) + 1;

            new File(workingDir + "/tmp").mkdir();
            List<RecordWriter> bucketWriters = new ObjectArrayList<RecordWriter>(numBuckets);
            for (int i = 0; i < numBuckets; i++) {
                bucketWriters.add(new RecordWriter(workingDir + "/tmp/bucket" + i, arguments.chunkSizePerWriter));
            }
            RecordWriter allWriter = new RecordWriter(args().outputFile);
            Random rand = new XoRoShiRo128PlusRandom();

            //set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "read";
            pgRead.expectedUpdates = totalRecords;
            pgRead.displayFreeMemory = true;
            pgRead.start();

            //fill buckets randomly
            System.out.println("Filling " + numBuckets + " temp buckets randomly");
            for (String filename : args().inputFiles) {
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
            System.out.printf("There are %d buckets to shuffle\n",numBuckets);
            //iterate over buckets
            ProgressLogger pgTempBucket = new ProgressLogger(LOG);
            pgTempBucket.itemsName = "buckets";
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



    @Override
    public RandomizerArguments createArguments() {
        return new RandomizerArguments();
    }



}
