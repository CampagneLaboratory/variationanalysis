package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.fastutil.ints.Int2BooleanAVLTreeMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2BooleanAVLTreeMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectArrayMap;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

/**
 * Combine a raw SBI with gold standard annotations to set isMutated flag. TSV format is chromosome\tposition\t[toBases]
 * where toBases is optional and contain the mutated bases (SNP or indel genotype of the mutation).
 * Created by fac2003 on 11/22/16.
 */
public class CombineWithGoldStandard extends AbstractTool<CombineWithGoldStandardArguments> {
    @Override
    public CombineWithGoldStandardArguments createArguments() {
        return new CombineWithGoldStandardArguments();
    }

    static private Logger LOG = LoggerFactory.getLogger(CombineWithGoldStandard.class);

    public static void main(String[] args) {
        CombineWithGoldStandard tool = new CombineWithGoldStandard();
        tool.parseArguments(args, "CombineWithGoldStandard", tool.createArguments());
        tool.execute();
    }

    @Override
    public void execute() {
        loadAnnotations(args().annotationFilename);
        RecordWriter outputWriters = null;
        try {
            outputWriters = new RecordWriter(args().outputFilename);
        } catch (IOException e) {
            throw new RuntimeException("Unable to create output file+" + args().outputFilename, e);
        }
        Random rand = new XorShift1024StarRandom();
        try (RecordReader reader = new RecordReader(args().sbiFilename)) {


//set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "records";
            pgRead.expectedUpdates = Math.min(reader.getTotalRecords(), reader.getTotalRecords());
            pgRead.displayFreeMemory = true;
            pgRead.start();
            long numWritten = 0;
            long numMutatedWritten = 0;
            for (BaseInformationRecords.BaseInformation record : reader) {
                record = annotate(record);
                double choice = rand.nextDouble();
                if (choice < args().samplingFraction || record.getMutated()) {
                    outputWriters.writeRecord(record);
                    numWritten += 1;
                    if (record.getMutated()) {
                        numMutatedWritten++;
                       /* System.out.printf("%s:%d\n\tnormal: %s\ttumor:  %s %n", record.getReferenceId(), record.getPosition(),
                                record.getSamples(0).getFormattedCounts(),
                                record.getSamples(1).getFormattedCounts());*/
                    }
                }
                pgRead.lightUpdate();
            }
            pgRead.stop();
            System.out.printf("Wrote %d records (of which %d are mutated).", numWritten, numMutatedWritten);
        } catch (IOException e) {
            throw new RuntimeException("Unable to load or write files. Check command line arguments.",e);
        } finally {
            IOUtils.closeQuietly(outputWriters);
        }
    }

    class Mutation {
        String from;
        String to;
        int position;
    }

    Object2ObjectArrayMap<String, Int2ObjectAVLTreeMap<Annotation>> annotations = new Object2ObjectArrayMap<>();

    class Annotation {
        int position;
        String toBase;
    }

    private void loadAnnotations(String annotationFilename) {
        try {
            LineIterator lines = new LineIterator(new FastBufferedReader(new FileReader(annotationFilename)));
            for (MutableString line : lines.allLines()) {
               if (line.startsWith("#")) continue;
                String tokens[] = line.toString().split("\t");
                final String chromosome = tokens[0];
                int position = Integer.parseInt(tokens[1]);
                // convert to zero-based position used by goby/variationanalysis:
                position -= 1;
                Int2ObjectAVLTreeMap<Annotation> perPositionMap = annotations.getOrDefault(chromosome, new Int2ObjectAVLTreeMap<>());                Annotation a = new Annotation();
                a.position = position;
                a.toBase = tokens.length >= 3 ? tokens[2] : null;
                perPositionMap.put(position, a);
                annotations.put(chromosome, perPositionMap);
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to find annotation filename:" + annotationFilename, e);

        }
    }

    private BaseInformationRecords.BaseInformation annotate(BaseInformationRecords.BaseInformation record) {
        if (!isAnnotated(record.getReferenceId(), record.getPosition())) {
            return record;
        }
        final BaseInformationRecords.BaseInformation.Builder builder = record.toBuilder();
        Annotation a=getAnnotation(record.getReferenceId(), record.getPosition());
        if (a.toBase!=null) {
            builder.setMutatedBase(a.toBase);
        }
        return builder.setMutated(true).build();
    }

    private boolean isAnnotated(String referenceId, int position) {
        return annotations.getOrDefault(referenceId, new Int2ObjectAVLTreeMap<Annotation>()).getOrDefault(position, null)!=null;
    }
    private Annotation getAnnotation(String referenceId, int position) {
        return annotations.getOrDefault(referenceId, new Int2ObjectAVLTreeMap<Annotation>()).get(position);
    }
}
