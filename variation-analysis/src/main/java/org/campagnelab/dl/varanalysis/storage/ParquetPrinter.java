package org.campagnelab.dl.varanalysis.storage;


import com.codahale.metrics.MetricRegistryListener;
import com.google.protobuf.TextFormat;
import org.apache.parquet.hadoop.ParquetReader;
import org.campagnelab.dl.varanalysis.intermediaries.Mutator;
import org.campagnelab.dl.varanalysis.intermediaries.Randomizer;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.eclipse.jetty.util.IO;

import java.io.IOException;
import java.util.List;


/**
 * Currently holds the main method. Jar takes two arguments
 * <p>
 * java -jar var-analysis.jar process /path/to/genotypes.parquet
 * creates mutated and randomized parquet file, then prints the latter
 * <p>
 * java -jar var-analysis.jar print /path/to/genotypes.parquet
 * prints the parquet file as is
 * <p>
 * Also, this jar should be a resource for Goby to output variations as a parquet file using its AvroVariationOutputFormat
 * class.
 * Created by rct66 on 5/17/16.
 *
 * @author rct66
 */
public class ParquetPrinter {

    String path;
    boolean focusPrint = false;
    private int refIndex;
    private int position;

    public static void main(String[] args) throws IOException {
        if (args.length<1) {
            System.err.println("usage: printer <parquet-file> [focus-ref-index focus-position]");
            System.exit(1);
        }
        ParquetPrinter parquetPrinter = new ParquetPrinter(args[0]);
        if (args.length >= 3) {
            // will only print the record(s) matching a specific position:
            parquetPrinter.setFocusOnPosition(Integer.parseInt(args[1]), Integer.parseInt(args[2]));
            System.out.println("Scanning for ");
        }
        parquetPrinter.print();

    }

    private void setFocusOnPosition(int refIndex, int position) {
        focusPrint = true;
        this.refIndex = refIndex;
        this.position = position;
    }


    public ParquetPrinter(String path) {
        this.path = path;
    }

    private void recordPrinter(BaseInformationRecords.BaseInformation base) throws IOException {
        TextFormat.print(base, System.out);
    }

    public void print() {
        try {
            RecordReader reader = new RecordReader(path);
            for (BaseInformationRecords.BaseInformation base : reader) {
                if (!focusPrint ||
                        (base.getReferenceIndex() == refIndex && base.getPosition() == position)) {
                    recordPrinter(base);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


}