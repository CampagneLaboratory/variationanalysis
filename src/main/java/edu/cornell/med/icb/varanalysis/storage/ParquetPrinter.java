package edu.cornell.med.icb.varanalysis.storage;



import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;

import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 */
public class ParquetPrinter {

    String path;

    public static void main(String args[]){
        ParquetPrinter printer = new ParquetPrinter(args[0]);
        printer.print();
    }

    public ParquetPrinter(String path){
        this.path = path;
    }

    private void recordPrinter(PosRecord prec){
        System.out.println("\n" + prec.getPosition() + "\t" + prec.getRefIdx() + prec.getMutated());
        List<SampleRecord> srecs = prec.getSamples();
        for (SampleRecord srec : srecs) {
            System.out.println(srec.getCounts().toString());
        }
    }

    public void print() {
        AvroVariationParquetReader avpReader = new AvroVariationParquetReader(path);
        while (true) {
            PosRecord pos = avpReader.read();
            if (pos == null) {
                break;
            }
            recordPrinter(pos);
        }
    }
}