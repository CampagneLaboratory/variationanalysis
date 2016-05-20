package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;

import java.io.File;

/**
 * Created by rct66 on 5/19/16.
 *
 * Describes an intermediary step in the sequence to prepare a parquet file for
 * deeplearning training
 *
 * All parquet files should uphold the same schema
 *
 */
public abstract class Intermediary {

    //pagesize may need to be bigger, but performance has been decent.
    static int blockSize = 256 * 1024 * 1024;
    static int pageSize = 64 * 1024;

    public void execute(String inPath, String outPath){
        execute(inPath,outPath,blockSize,pageSize);
    }

    public abstract void execute(String inPath, String outPath, int blockSize, int pageSize);


}