package org.campagnelab.dl.somatic.tools;


import com.google.protobuf.TextFormat;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;

/**
 * A tool to print an sbi to text.
 * <p>
 * Created by Fabien Campagne on 7/23/17.
 *
 * @author Fabien Campagne
 */
public class Print extends AbstractTool<PrintArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Print.class);

    public static void main(String[] args) {

        Print tool = new Print();
        tool.parseArguments(args, "Print", tool.createArguments());
        tool.execute();
    }


    @Override
    public void execute() {

        try {
            long totalRecords = 0;

            RecordReader source = new RecordReader(args().inputFile);
            totalRecords = source.getTotalRecords();


            source = new RecordReader(args().inputFile);
            for (BaseInformationRecords.BaseInformation inputRecord : source) {

                TextFormat.print(inputRecord, System.out);

            }
            source.close();


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
        if (outputComplexityLevel >= 2) {

            System.err.printf("output complexity level %d is not supported at this time.%n", outputComplexityLevel);
            System.exit(1);

        }
        return null;
    }


    @Override
    public PrintArguments createArguments() {
        return new PrintArguments();
    }


}
