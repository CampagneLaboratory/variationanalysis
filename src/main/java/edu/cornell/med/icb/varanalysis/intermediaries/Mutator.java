package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;

/**
 * Created by rct66 on 5/18/16.
 *
 * the mutator object iterates over a parquet file and creates am additional mutated copy of every record.
 */
public class Mutator extends Intermediary{

    //min fraction of bases mutated at a record
    int deltaSmall;
    //max fraction of bases mutated at a record
    int deltaBig;


    @Override
    void Process(AvroVariationParquetReader reader, AvroVariationParquetWriter writer) {

    }
}
