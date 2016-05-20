package edu.cornell.med.icb.varanalysis.storage;



import edu.cornell.med.icb.varanalysis.format.PosRecord;
import edu.cornell.med.icb.varanalysis.format.SampleRecord;
import edu.cornell.med.icb.varanalysis.intermediaries.Intermediary;
import edu.cornell.med.icb.varanalysis.intermediaries.Mutator;
import edu.cornell.med.icb.varanalysis.intermediaries.Randomizer;

import java.util.List;


/**
 * Created by rct66 on 5/17/16.
 *
 * Currently holds the main method. Jar takes two arguments
 *
 *  java -jar var-analysis.jar process /path/to/genotypes.parquet
 *  creates mutated and randomized parquet file, then prints the latter
 *
 *  java -jar var-analysis.jar print /path/to/genotypes.parquet
 *  prints the parquet file as is
 *
 * Also, this jar should be a resource for Goby to output variations as a parquet file using its AvroVariationOutputFormat
 * class.
 *
 */
public class ParquetPrinter {

    String path;

    public static void main(String[] args){
        if (args[0].equals("process")){
            main2(args[1]);
        } else {
            ParquetPrinter printer = new ParquetPrinter(args[1]);
            printer.print();
        }
    }

    public static void main2(String startPath){

        //mutate
        Mutator mut = new Mutator();
        String mutPath = startPath.substring(0, startPath.length() - 8) + "_mutated.parquet";
        mut.execute(startPath,mutPath);
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