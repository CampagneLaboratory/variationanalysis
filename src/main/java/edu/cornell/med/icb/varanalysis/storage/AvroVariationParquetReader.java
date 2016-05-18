package edu.cornell.med.icb.varanalysis.storage;



import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.avro.AvroParquetReader;
import org.apache.parquet.hadoop.ParquetReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 *
 * This is the avp reader class. All changes to the posrecord or samplerecord classes (or record.avsc schema)
 * must be reflected in this reader.
 * This reader serves as a consistent interface with those changing storage models.
 *
 *
 */



public class AvroVariationParquetReader {

    String path;
    ParquetReader<GenericRecord> pqReader;

    public AvroVariationParquetReader(String path){
        this.path = path;
        AvroParquetReader.Builder<GenericRecord> pqBuilder = AvroParquetReader.builder(new Path(path));
        try {
            pqReader = pqBuilder.build();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    public void close(){
        try {
            pqReader.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public PosRecord read(){
        try {
            //hacky workaround to classcastexception raised when casting GenericRecords as PosRecords
            //TODO: when there is time, look for classpath fix

            GenericRecord pos = pqReader.read();
            if (pos == null){
                return null;
            }
            PosRecord specPos = new PosRecord();
            specPos.setPosition((Integer)pos.get("position"));
            specPos.setRefIdx((Integer)pos.get("refIdx"));
            GenericData.Array<GenericRecord> samples = (GenericData.Array<GenericRecord>)pos.get("samples");
            List<SampleRecord> sampleList= new ArrayList<SampleRecord>();
            for (GenericRecord sampleG : samples){
                SampleRecord sample = new SampleRecord();
                List<Integer> countsList = new ArrayList<Integer>();
                GenericData.Array<Integer> counts = (GenericData.Array<Integer>)sampleG.get("counts");
                for (Integer i : counts){
                    countsList.add(i);
                }
                sample.setCounts(countsList);
                sampleList.add(sample);
            }
            specPos.setSamples(sampleList);
            return specPos;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}