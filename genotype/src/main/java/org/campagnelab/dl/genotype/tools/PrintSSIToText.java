package org.campagnelab.dl.genotype.tools;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.storage.SegmentReader;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * A tool to print an ssi to text.
 *
 * @author manuele
 */
public class PrintSSIToText extends AbstractTool<PrintSSIToTextArguments> {

    static private final Logger LOG = LoggerFactory.getLogger(PrintSSIToText.class);

    public static void main(String[] args) {

        PrintSSIToText tool = new PrintSSIToText();
        tool.parseArguments(args, "PrintSSIToText", tool.createArguments());
        tool.execute();
    }

    @Override
    public PrintSSIToTextArguments createArguments() {
        return new PrintSSIToTextArguments();
    }

    @Override
    public void execute() {
        try {
            long totalRecords = 0;

            SegmentReader source = new SegmentReader(args().inputFile);
            totalRecords = source.getTotalRecords();

            source = new SegmentReader(args().inputFile);

            for (SegmentInformationRecords.SegmentInformation segment : source) {
                if (args().removeFeatures || args().removeLabels) {
                    SegmentInformationRecords.SegmentInformation.Builder segBuilder = segment.toBuilder();
                    int sampleIndex = 0;
                    for (SegmentInformationRecords.Sample sample : segBuilder.getSampleList()) {
                        SegmentInformationRecords.Sample.Builder sampleBuilder = sample.toBuilder();
                        int baseIndex = 0;

                        for (SegmentInformationRecords.Base base : sample.getBaseList()) {
                            SegmentInformationRecords.Base.Builder baseBuilder = base.toBuilder();
                            if (args().removeFeatures) {
                                baseBuilder.clearFeatures();
                            }
                            if (args().removeLabels) {
                                baseBuilder.clearLabels();
                            }
                            sampleBuilder.setBase(baseIndex, baseBuilder);
                            baseIndex++;
                        }
                        segBuilder.setSample(sampleIndex, sampleBuilder);
                    }
                    segment=segBuilder.build();
                }
                TextFormat.print(segment, System.out);

            }
            source.close();
            System.out.println("Total records: " + totalRecords);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
