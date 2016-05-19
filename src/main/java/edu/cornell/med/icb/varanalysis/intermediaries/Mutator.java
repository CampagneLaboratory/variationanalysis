package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;

/**
 * Created by rct66 on 5/18/16.
 *
 * the mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 */
public class Mutator{

    //min fraction of bases mutated at a record
    double deltaSmall = 0.05;
    //max fraction of bases mutated at a record
    double deltaBig = 0.5;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;



    void Process(String in, String out, int blockSize, int pageSize) {
        AvroVariationParquetReader reader = new AvroVariationParquetReader(in);
        AvroVariationParquetWriter writer = new AvroVariationParquetWriter(out,blockSize,pageSize);
    }

    int[] Mutate(int[] forward, int[] backward)
}
