package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;

/**
 * Created by rct66 on 5/18/16.
 *
 * the randomizer object iterates over a parquet file and randomizes the order of records in batches
 */
public class Randomizer extends Intermediary{


    @Override
    void Process(AvroVariationParquetReader reader, AvroVariationParquetWriter writer) {

    }
}
