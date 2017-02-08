package org.campagnelab.dl.somatic.storage;


import com.google.protobuf.TextFormat;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Random;
import java.util.Set;


/**
 *
 * Simple downsampler to create an sbi with some proportion of input records.
 * Created by rct66 on 2/2/16.
 *
 * @author rct66
 */
public class DownSampler {

    private static final Logger LOG = LoggerFactory.getLogger(DownSampler.class);

    String path;
    static int actualCount = 0;
    private boolean makeDebug = true;
    private RecordWriter makeDebugWriter;
    Random rand = new XoRoShiRo128PlusRandom();
    Float rate;

    public static void main(String[] args) throws IOException {
        if (args.length<1) {
            System.err.println("usage: printer <parquet-file> [focus-ref-index focus-position]");
            System.exit(1);
        }
        DownSampler protobufPrinter = new DownSampler(args[0],args[1]);

        protobufPrinter.downs();
        System.out.println("actual count: " + actualCount);

    }




    public DownSampler(String path, String rate) throws IOException {
        this.path = path;
        this.rate = Float.parseFloat(rate);
        if (makeDebug){
            String debugPath = FilenameUtils.getFullPath(path) + FilenameUtils.getBaseName(path) + "_" + rate.toString() + ".sbi";
            makeDebugWriter = new RecordWriter(debugPath);
        }
    }


    public void downs() {
        ProgressLogger recordLogger = new ProgressLogger(LOG);
        try {
            RecordReader reader = new RecordReader(path);
            System.out.println(reader.getTotalRecords());
            recordLogger.expectedUpdates = reader.getTotalRecords();
            for (BaseInformationRecords.BaseInformation base : reader) {
                if (rand.nextFloat() <= rate){
                    makeDebugWriter.writeRecord(base);
                    actualCount++;
                }
                recordLogger.lightUpdate();

            }
            makeDebugWriter.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


}