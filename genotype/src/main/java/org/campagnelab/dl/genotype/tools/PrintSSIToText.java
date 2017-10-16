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

            for (SegmentInformationRecords.SegmentInformation segment : source)
                TextFormat.print(segment, System.out);
            source.close();
            System.out.println("Total records: " + totalRecords);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
