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
    int blockSize = 256 * 1024 * 1024;
    int pageSize = 64 * 1024;

    public void execute(String in, String out){
        AvroVariationParquetReader reader = new AvroVariationParquetReader(in);
        AvroVariationParquetWriter writer = new AvroVariationParquetWriter(out,blockSize,pageSize);
        Process(reader,writer);
    }

    abstract void Process(AvroVariationParquetReader reader, AvroVariationParquetWriter writer);


}
