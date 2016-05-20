package edu.cornell.med.icb.varanalysis.storage;



import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import edu.cornell.med.icb.varanalysis.intermediaries.Intermediary;
import edu.cornell.med.icb.varanalysis.intermediaries.Mutator;

import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 */
public class ParquetPrinter {

    String path;

    public static void main(String args[]){
        Mutator mut = new Mutator();
        String mutPath = args[0]+"_mutated";
        mut.process(args[0],args[0]+"_mutated", Intermediary.blockSize, Intermediary.pageSize);
        ParquetPrinter printer = new ParquetPrinter(mutPath);
        printer.print();
    }

    public ParquetPrinter(String path){
        this.path = path;
    }

    private void recordPrinter(PosRecord prec){
        System.out.println("\n" + prec.getPosition() + "\t" + prec.getRefIdx() + "\t" + prec.getMutated());
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