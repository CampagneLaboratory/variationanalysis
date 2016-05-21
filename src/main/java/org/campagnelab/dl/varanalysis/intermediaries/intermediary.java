package org.campagnelab.dl.varanalysis.intermediaries;

/**
 *
 *
 * Describes an intermediary step in the sequence to prepare a parquet file for
 * deeplearning training
 *
 * All parquet files should uphold the same schema
 * Created by rct66 on 5/19/16.
 * @author rct66
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