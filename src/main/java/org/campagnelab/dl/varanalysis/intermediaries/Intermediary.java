package org.campagnelab.dl.varanalysis.intermediaries;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

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

    public void executeOver(String inPath, String outPath) throws IOException {
        java.nio.file.Path toDelete = Paths.get(inPath);
        Files.deleteIfExists(toDelete);
        execute(inPath,outPath,blockSize,pageSize);
    }

    public void executeOver(String inPath, String outPath, int blockSize, int pageSize) throws IOException {
        java.nio.file.Path toDelete = Paths.get(inPath);
        Files.deleteIfExists(toDelete);
        execute(inPath,outPath,blockSize,pageSize);
    }

    public void execute(String inPath, String outPath) throws IOException {
        execute(inPath,outPath,blockSize,pageSize);
    }

    public abstract void execute (String inPath, String outPath, int blockSize, int pageSize) throws IOException ;


}