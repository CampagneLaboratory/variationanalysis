package edu.cornell.med.icb.varanalysis.intermediaries;

import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetReader;
import edu.cornell.med.icb.varanalysis.storage.AvroVariationParquetWriter;

import java.util.List;

/**
 * Created by rct66 on 5/18/16.
 *
 * the mutator object iterates over a parquet file and creates an additional mutated copy of every record.
 * Mutations
 *
 */
public class Mutator{
    //delta will be halved in homozygous cases (to account for twice the reads at a base)
    //min fraction of bases mutated at a record
    double deltaSmall = 0.1;
    //max fraction of bases mutated at a record
    double deltaBig = 1.0;
    //minimum proportion of counts to presume allele
    double zygHeuristic = 0.1;



    void process(String in, String out, int blockSize, int pageSize) {
        AvroVariationParquetReader reader = new AvroVariationParquetReader(in);
        AvroVariationParquetWriter writer = new AvroVariationParquetWriter(out,blockSize,pageSize);
        PosRecord pos;
        while(true) {
            //read and write unmutated record
            pos = reader.read();
            if (pos == null) break;
            writer.writeRecord(pos);
            //mutate record and write it again
            pos.setMutated(true);
            List<SampleRecord> sampleList = pos.getSamples();
            SampleRecord mutating = sampleList.get(sampleList.size() - 1);
            mutating.setCounts(mutate(mutating.getCounts()));
            pos.setSamples(sampleList);
            //record classes do not deep copy
            writer.writeRecord(pos);
        }

    }

    //backward appended to forward in input and output
    //10 fields: ATCGOATCGO
    //List instead of array because of avro code generation...
    private List<Integer> mutate(List<Integer> counts){
        return null;
    }

}
