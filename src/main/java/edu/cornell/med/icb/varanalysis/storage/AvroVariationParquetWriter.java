package edu.cornell.med.icb.varanalysis.storage;

import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import org.apache.avro.Schema;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.avro.AvroParquetWriter;
import org.apache.parquet.avro.AvroSchemaConverter;
import org.apache.parquet.avro.AvroWriteSupport;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.MessageType;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by rct66 on 5/17/16.
 *
 *
 *
 * This is the avp writer class. All changes to the posrecord or samplerecord classes (or record.avsc schema)
 * must be reflected in this writer.
 * This reader serves as a consistent interface with those changing storage models.
 *
 *
 *
 */
public class AvroVariationParquetWriter {
    ParquetWriter parquetWriter;

    public AvroVariationParquetWriter(String path, int blockSize, int pageSize){
        // load your Avro schema

        Schema avroSchema = null;
        try {
            avroSchema = new Schema.Parser().parse(this.getClass().getClassLoader().getResourceAsStream("record.avsc"));

            // generate the corresponding Parquet schema
            MessageType parquetSchema = new AvroSchemaConverter().convert(avroSchema);

            // create a WriteSupport object to serialize your Avro objects
            AvroWriteSupport writeSupport = new AvroWriteSupport(parquetSchema, avroSchema);

            // choose compression scheme
            //compressionCodecName = CompressionCodecName.SNAPPY;


            Path outputPath = new Path(path);

            // the ParquetWriter object that will consume Avro GenericRecords
            parquetWriter = new AvroParquetWriter(outputPath,
                    avroSchema, CompressionCodecName.UNCOMPRESSED, blockSize, pageSize);


        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    public void writeRecord(PosRecord pos){
        try {
            parquetWriter.write(pos);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void writeRecord(List<int[][]> samplesCounts, int referenceIndex, int position, boolean mutated) {

        PosRecord pos = new PosRecord();
        pos.setPosition(position);
        pos.setRefIdx(referenceIndex);
        pos.setMutated(mutated);
        List<SampleRecord> samples = new ArrayList<SampleRecord>();
        for (int[][] countsArr : samplesCounts){
            SampleRecord rec = new SampleRecord();
            List<Integer> countsList = new ArrayList<Integer>();
            for(int j=0;j<2;j++)
                for(int k=0;k<5;k++)
                    countsList.add(countsArr[j][k]);
            rec.setCounts(countsList);
            samples.add(rec);
        }
        pos.setSamples(samples);
        writeRecord(pos);

    }

    public void close() {
        try {
            parquetWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
