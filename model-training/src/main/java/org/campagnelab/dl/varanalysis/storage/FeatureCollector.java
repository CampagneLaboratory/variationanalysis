package org.campagnelab.dl.varanalysis.storage;


import com.fasterxml.jackson.core.util.BufferRecycler;
import com.google.protobuf.TextFormat;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV18;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Set;


/**
 * Currently holds the main method. Jar takes two arguments
 * <p>
 * java -jar var-analysis.jar process /recordsPath/to/genotypes.parquet
 * creates mutated and randomized parquet file, then prints the latter
 * <p>
 * java -jar var-analysis.jar print /recordsPath/to/genotypes.parquet
 * prints the parquet file as is
 * <p>
 * Also, this jar should be a resource for Goby to output variations as a parquet file using its AvroVariationOutputFormat
 * class.
 * Created by rct66 on 5/17/16.
 *
 * @author rct66
 */
public class FeatureCollector {
    private static FeatureMapper mapper = new FeatureMapperV18();
    String recordsPath;
    BufferedReader positionsReader;
    private Set<Integer> posSet = new IntOpenHashSet();

    public static void main(String[] args) throws IOException {
        if (args.length<1) {
            System.err.println("usage: printer <positions-file> <parquet-file>");
            System.exit(1);
        }
        FeatureCollector featureCollector = new FeatureCollector(args[0], args[1]);
        featureCollector.getPositions();
        featureCollector.print();



    }

    public FeatureCollector(String positionsPath, String recordsPath) throws FileNotFoundException {
        this.recordsPath = recordsPath;
        this.positionsReader =  new BufferedReader((new FileReader(positionsPath)));
    }

    private void recordPrinter(BaseInformationRecords.BaseInformation base) throws IOException {
        //TODO
    }

    private void getPositions() throws IOException{
        String line;
        while ((line = positionsReader.readLine()) != null) {
            posSet.add(Integer.parseInt(line));
        }
    }

    public void print() {
        try {
            RecordReader reader = new RecordReader(recordsPath);
            for (BaseInformationRecords.BaseInformation base : reader) {
                if (posSet.contains(base.getPosition())) {
                    recordPrinter(base);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


}