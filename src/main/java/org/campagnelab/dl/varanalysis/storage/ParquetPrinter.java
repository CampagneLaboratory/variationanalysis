package org.campagnelab.dl.varanalysis.storage;




import com.codahale.metrics.MetricRegistryListener;
import com.google.protobuf.TextFormat;
import org.apache.parquet.hadoop.ParquetReader;
import org.campagnelab.dl.varanalysis.intermediaries.Randomizer;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.IOException;
import java.util.List;


/**
 *
 *
 * Currently holds the main method. Jar takes two arguments
 *
 *  java -jar var-analysis.jar process /path/to/genotypes.parquet
 *  creates mutated and randomized parquet file, then prints the latter
 *
 *  java -jar var-analysis.jar print /path/to/genotypes.parquet
 *  prints the parquet file as is
 *
 * Also, this jar should be a resource for Goby to output variations as a parquet file using its AvroVariationOutputFormat
 * class.
 * Created by rct66 on 5/17/16.
 * @author rct66
 */
public class ParquetPrinter {

    String path;

    public static void main(String[] args){
        ParquetPrinter printer = new ParquetPrinter(args[0]);
        printer.print();
    }


    public ParquetPrinter(String path){
        this.path = path;
    }

    private void recordPrinter(BaseInformationRecords.BaseInformation base) throws IOException{
        TextFormat.print(base, System.out);
    }

    public void print() {
        try {
            RecordReader reader = new RecordReader(path);
            for (BaseInformationRecords.BaseInformation base : reader){
                recordPrinter(base);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}