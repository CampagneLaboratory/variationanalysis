package org.campagnelab.dl.varanalysis.intermediaries;


import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetReader;
import org.campagnelab.dl.varanalysis.storage.AvroVariationParquetWriter;
import it.unimi.dsi.util.XorShift128PlusRandom;
import org.apache.avro.generic.GenericRecord;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Created by rct66 on 5/18/16.
 *
 * the randomizer object iterates over a parquet file and randomizes the order of records in batches
 */
public class Randomizer extends Intermediary{

    public void execute(String inPath, String outPath, int blockSize, int pageSize) {
        AvroVariationParquetReader reader = new AvroVariationParquetReader(inPath);
        AvroVariationParquetWriter writer = new AvroVariationParquetWriter(outPath, blockSize, pageSize);
        GenericRecord gen;
        List<GenericRecord> genList = new ArrayList<GenericRecord>();
        while (true){
            gen = reader.read();
            if (gen == null)
                break;
            genList.add(gen);
        }
        reader.close();

        Random rand = new XorShift128PlusRandom();
        Collections.shuffle(genList,rand);

        for (GenericRecord outgen : genList){
            writer.writeRecord(outgen);
        }
        writer.close();
    }
}
