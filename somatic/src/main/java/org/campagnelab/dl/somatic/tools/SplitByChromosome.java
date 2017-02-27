package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import play.mvc.WebSocket;

import java.io.IOException;
import java.util.Map;
import java.util.Set;

/**
 * Split a BSI file into several parts. Useful for creating training/validation/test splits of a large dataset.split
 * Created by fac2003 on 9/2/16.
 */
public class SplitByChromosome extends AbstractTool<SplitByChromosomeArguments> {


    static private Logger LOG = LoggerFactory.getLogger(SplitByChromosome.class);

    public static void main(String[] args) {

        SplitByChromosome tool = new SplitByChromosome();
        tool.parseArguments(args, "Split", tool.createArguments());
        tool.execute();
    }
    @Override
    public void execute()  {

        Set<String> testIDs = new ObjectArraySet<>(args().testChromosomes);
        Set<String> valIDs = new ObjectArraySet<>(args().valChromosomes);
        Map<String,Integer> testCounts = new Object2IntArrayMap<>(testIDs.size());
        Map<String,Integer> valCounts = new Object2IntArrayMap<>(valIDs.size());
        Map<String,Integer> trainCounts = new Object2IntArrayMap<>(10);



        try (RecordReader reader = new RecordReader(args().inputFile)) {

            RecordWriter trainWriter = new RecordWriter(args().outputFile + "train");
            RecordWriter valWriter = new RecordWriter(args().outputFile + "validation");
            RecordWriter testWriter = new RecordWriter(args().outputFile + "test");

            //set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "records";
            pgRead.expectedUpdates = reader.getTotalRecords();
            pgRead.displayFreeMemory = true;
            pgRead.start();
            long numWritten = 0;

            for (BaseInformationRecords.BaseInformation record : reader) {
                String refID = record.getReferenceId();
                //increment count of chrom for print output


                //write record to appropriate writer
                if (testIDs.contains(refID)){
                    testWriter.writeRecord(record);
                    int count = testCounts.computeIfAbsent(refID, k -> 0);
                    testCounts.put(refID,++count);
                } else if (valIDs.contains(refID)){
                    valWriter.writeRecord(record);
                    int count = valCounts.computeIfAbsent(refID, k -> 0);
                    valCounts.put(refID,++count);
                } else {
                    trainWriter.writeRecord(record);
                    int count = trainCounts.computeIfAbsent(refID, k -> 0);
                    trainCounts.put(refID,++count);
                }
                pgRead.lightUpdate();
                numWritten += 1;

            }
            pgRead.stop();
            trainWriter.close();
            testWriter.close();
            valWriter.close();

            int sumTrain = trainCounts.values().stream().mapToInt(Integer::intValue).sum();
            int sumVal = valCounts.values().stream().mapToInt(Integer::intValue).sum();
            int sumTest = testCounts.values().stream().mapToInt(Integer::intValue).sum();
            float fractionTrain = (float)sumTrain/reader.numRecords();
            float fractionVal = (float)sumVal/reader.numRecords();
            float fractionTest = (float)sumTest/reader.numRecords();

            System.out.println("train counts = " + sumTrain + "," + fractionTrain + ": " + trainCounts);
            System.out.println("validation counts = " + sumVal + "," + fractionVal + ": " + valCounts);
            System.out.println("test counts = " + sumTest + "," + fractionTest + ": " + testCounts);

            reader.close();

        } catch (IOException e) {
            System.err.println("Unable to load or write files. Check command line arguments.");
        }
    }

   

    @Override
    public SplitByChromosomeArguments createArguments() {
            return new SplitByChromosomeArguments();
    }


}
