package org.campagnelab.dl.somatic.tools;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.apache.commons.io.FileUtils;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * The simplify tool removes fields from the sbi file to produce simpler messages. Use to create
 * a progressive curiculum for DL networks to learn an sbi to answer mapping directly.
 * <p>
 * Created by Fabien Campagne on 7/23/17.
 *
 * @author Fabien Campagne
 */
public class Simplify extends AbstractTool<SimplifyArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Simplify.class);

    public static void main(String[] args) {

        Simplify tool = new Simplify();
        tool.parseArguments(args, "Simplify", tool.createArguments());
        tool.execute();
    }


    @Override
    public void execute() {
        String workingDir = new File(args().outputFile).getParent();
        if (workingDir == null) {
            workingDir = ".";
        }
        try {
            long totalRecords = 0;

            RecordReader source = new RecordReader(args().inputFile);
            totalRecords = source.getTotalRecords();

            System.out.println("Simplifying input to level "+args().outputComplexityLevel);
            RecordWriter allWriter = new RecordWriter(args().outputFile);

            //set up logger
            ProgressLogger pgRead = new ProgressLogger(LOG);
            pgRead.itemsName = "records";
            pgRead.expectedUpdates = totalRecords;
            pgRead.displayFreeMemory = true;
            pgRead.start();


            source = new RecordReader(args().inputFile);
            for (BaseInformationRecords.BaseInformation inputRecord : source) {

                allWriter.writeRecord(simplify(args().outputComplexityLevel, inputRecord));
                pgRead.lightUpdate();
            }
            source.close();
            allWriter.close();
            System.out.println("Simplified file written to "+args().outputFile+".sbi");
            pgRead.stop();

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private BaseInformationRecords.BaseInformation simplify(int outputComplexityLevel, BaseInformationRecords.BaseInformation inputRecord) {
        BaseInformationRecords.BaseInformation.Builder builder = inputRecord.toBuilder();
        if (outputComplexityLevel == 0) {
            // level 0 removes all samples:
            builder.clearSamples();
        } else {
            int sampleIndex = 0;
            for (BaseInformationRecords.SampleInfo sample : builder.getSamplesList()) {
                builder.setSamples(sampleIndex++, simplify(outputComplexityLevel, sample));
            }
        }

        return builder.build();
    }

    private BaseInformationRecords.SampleInfo simplify(int outputComplexityLevel, BaseInformationRecords.SampleInfo sample) {
        BaseInformationRecords.SampleInfo.Builder builder = sample.toBuilder();
        if (outputComplexityLevel == 1) {
            // level 1 removes all counts:
            builder.clearCounts();
        } else {
            int countIndex = 0;
            for (BaseInformationRecords.CountInfo count : builder.getCountsList()) {
                builder.setCounts(countIndex++, simplify(outputComplexityLevel, count));
            }
        }
        return builder.build();
    }

    private BaseInformationRecords.CountInfo simplify(int outputComplexityLevel, BaseInformationRecords.CountInfo count) {

        BaseInformationRecords.CountInfo.Builder builder = count.toBuilder();
        if (outputComplexityLevel == 2) {
            // level 2 removes all frequencies:
            builder.clearDistancesToReadVariationsForwardStrand();
            builder.clearDistancesToReadVariationsReverseStrand();
            builder.clearDistanceToEndOfRead();
            builder.clearDistanceToStartOfRead();
            builder.clearInsertSizes();
            builder.clearQueryAlignedLengths();
            builder.clearTargetAlignedLengths();
            builder.clearQualityScoresForwardStrand();
            builder.clearQualityScoresReverseStrand();
            builder.clearNumVariationsInReads();
            builder.clearReadMappingQualityForwardStrand();
            builder.clearReadMappingQualityReverseStrand();
        } else {
            if (outputComplexityLevel >= 3) {

                System.err.printf("output complexity level %d is not supported at this time.%n", outputComplexityLevel);
                System.exit(1);

            }
        }
        return builder.build();

    }


    @Override
    public SimplifyArguments createArguments() {
        return new SimplifyArguments();
    }


}
