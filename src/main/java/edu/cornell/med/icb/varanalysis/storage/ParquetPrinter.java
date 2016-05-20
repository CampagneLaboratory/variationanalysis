package edu.cornell.med.icb.varanalysis.storage;



import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import edu.cornell.med.icb.varanalysis.intermediaries.Intermediary;
import edu.cornell.med.icb.varanalysis.intermediaries.Mutator;
import edu.cornell.med.icb.varanalysis.intermediaries.Randomizer;

import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 */
public class ParquetPrinter {

    String path;

    public static void main(String args[]){


        String startPath = args[0];

        //mutate
        Mutator mut = new Mutator();
        String mutPath = startPath.substring(0, startPath.length() - 8) + "_mutated.parquet";
        mut.execute(args[0],mutPath);
        System.out.println("mutated");
        //randomize
        Randomizer rndz = new Randomizer();
        String rndPath = mutPath.substring(0, startPath.length() - 8) + "_randomized.parquet";
        rndz.execute(mutPath,rndPath);
        System.out.println("randomized");

        ParquetPrinter printer = new ParquetPrinter(rndPath);
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