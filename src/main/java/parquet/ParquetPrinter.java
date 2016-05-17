package parquet;



import avro.PosRecord;
import avro.SampleRecord;
import org.apache.avro.Schema;
import org.apache.avro.generic.IndexedRecord;
import org.apache.hadoop.record.Index;
import org.apache.parquet.avro.AvroReadSupport;
import org.apache.parquet.avro.AvroSchemaConverter;
import org.apache.parquet.hadoop.ParquetReader;
import org.apache.parquet.schema.MessageType;

import java.io.File;
import java.io.IOException;
import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 */
public class ParquetPrinter {

    public void recordPrinter(PosRecord prec){
        System.out.println("\n" + prec.getPosition() + "\t" + prec.getRefIdx());
        List<SampleRecord> srecs = prec.getSamples();
        for (SampleRecord srec : srecs) {
            System.out.println(srec.getCounts().toString());
        }
    }

    public void pqFilePrinter(File file){
        Schema avroSchema = null;
        try {
            avroSchema = new Schema.Parser().parse(ClassLoader.getSystemResourceAsStream("record.avsc"));


            // generate the corresponding Parquet schema
            MessageType parquetSchema = new AvroSchemaConverter().convert(avroSchema);

            // create a WriteSupport object to serialize your Avro objects
            //AvroWriteSupport writeSupport = new AvroWriteSupport(parquetSchema, avroSchema);
            AvroReadSupport<IndexedRecord> readSupport = new AvroReadSupport<IndexedRecord>();

            // choose compression scheme
            //compressionCodecName = CompressionCodecName.SNAPPY;

            // set Parquet file block size and page size values
            //int blockSize = 256 * 1024 * 1024;
            //int pageSize = 64 * 1024;

            //Path outputPath = new Path(outputFilename);

            // the ParquetWriter object that will consume Avro GenericRecords
            parquetWriter = new ParquetWriter(outputPath,
                    writeSupport, compressionCodecName, blockSize, pageSize);

            for (GenericRecord record : SourceOfRecords) {
                parquetWriter.write(record);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
