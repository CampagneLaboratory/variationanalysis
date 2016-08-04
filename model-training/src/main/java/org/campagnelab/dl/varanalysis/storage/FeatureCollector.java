package org.campagnelab.dl.varanalysis.storage;


import it.unimi.dsi.fastutil.objects.ObjectOpenHashSet;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV18;
import org.campagnelab.dl.model.utils.mappers.FeatureNameMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.io.*;
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
    private static FeatureNameMapper mapper = new FeatureMapperV18();
    String recordsPath;
    BufferedReader positionsReader;
    BufferedWriter outputWriter;
    private Set<String> idSet = new ObjectOpenHashSet<String>();

    public static void main(String[] args) throws IOException {
        if (args.length<1) {
            System.err.println("usage: printer <positions-file> <parquet-file> <feature-output>");
            System.exit(1);
        }
        FeatureCollector featureCollector = new FeatureCollector(args[0], args[1], args[2]);
        featureCollector.getIds();
        featureCollector.execute();



    }

    public FeatureCollector(String positionsPath, String recordsPath, String outputPath) throws IOException {
        this.recordsPath = recordsPath;
        this.positionsReader =  new BufferedReader((new FileReader(positionsPath)));
        this.outputWriter = new BufferedWriter(new FileWriter(outputPath));
    }

    private void outputRecord(BaseInformationRecords.BaseInformation base) throws IOException {
        String out =
                base.getMutated()?"1":"0" + "\t"
                + base.getReferenceIndex() + ":" + base.getPosition() + "\t"
                + getSumCounts(base)
                + getFeatures(base) + "\n";
        outputWriter.append(out);
    }

    private String getFeatures(BaseInformationRecords.BaseInformation base) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < mapper.numberOfFeatures(); i++){
            sb.append("\t" + mapper.produceFeature(base,i));
        }
        return sb.toString();
    }

    private String getSumCounts(BaseInformationRecords.BaseInformation base) {
        int sum = 0;
        for (BaseInformationRecords.SampleInfo sample : base.getSamplesList()){
            for (BaseInformationRecords.CountInfo count : sample.getCountsList()){
                sum += (count.getGenotypeCountForwardStrand() + count.getGenotypeCountReverseStrand());
            }
        }
        return Integer.toString(sum);
    }

    private void getIds() throws IOException{
        String line;
        while ((line = positionsReader.readLine()) != null) {
            idSet.add(line);
        }
    }

    private void execute() throws IOException {
        outputWriter.write(getHeader());
        try {
            RecordReader reader = new RecordReader(recordsPath);
            for (BaseInformationRecords.BaseInformation base : reader) {
                if (idSet.contains(base.getReferenceIndex()+":"+base.getPosition())) {
                    outputRecord(base);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    private String getHeader() {
        StringBuffer sb = new StringBuffer();
        sb.append("group\tid");
        for (int i = 0; i < mapper.numberOfFeatures(); i++){
            sb.append("\t" + mapper.getFeatureName(i));
        }
        return sb.toString();
    }
}