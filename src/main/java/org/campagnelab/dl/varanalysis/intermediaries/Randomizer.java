package org.campagnelab.dl.varanalysis.intermediaries;


import org.apache.hadoop.fs.FileUtil;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.campagnelab.dl.varanalysis.storage.RecordWriter;
import it.unimi.dsi.util.XorShift128PlusRandom;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 *
 *
 * the randomizer object iterates over a parquet file and randomizes the order of records in batches
 *
 * Created by rct66 on 5/18/16.
 * @author rct66
 */
public class Randomizer extends Intermediary{

    public static void main (String[] args){

        //randomize
        Randomizer rndz = new Randomizer();
        String startPath = args[0];
        String rndPath = args[0].substring(0, startPath.length() - 8) + "_randomized.parquet";
        rndz.execute(startPath,rndPath, blockSize, pageSize);
        System.out.println("randomized");
    }

    public void execute(String inPath, String outPath, int blockSize, int pageSize) {
        try {
            RecordReader reader = new RecordReader(inPath);

            RecordWriter writer = new RecordWriter(outPath,blockSize,pageSize);

            List<BaseInformationRecords.BaseInformation> recList = new ArrayList<BaseInformationRecords.BaseInformation>();
            for (BaseInformationRecords.BaseInformation rec : reader) {
                recList.add(rec);
            }
            reader.close();

            Random rand = new XorShift128PlusRandom();
            Collections.shuffle(recList,rand);

            for ( BaseInformationRecords.BaseInformation rec : recList){
                writer.writeRecord(rec);
            }
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
