package org.campagnelab.dl.somatic.tools;


import com.google.protobuf.TextFormat;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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

                TextFormat.print(simplify(args().simplifyLevel, inputRecord), System.out);

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
            ObjectArrayList<BaseInformationRecords.CountInfo> toAdd=new ObjectArrayList<>();
            for (BaseInformationRecords.CountInfo count : builder.getCountsList()) {
                final BaseInformationRecords.CountInfo simplified = simplify(outputComplexityLevel, count);
                if (simplified!=null) {
                    toAdd.add(simplified);
                }
            }
            builder.clearCounts();
            builder.addAllCounts(toAdd);
        }
        return builder.build();
    }

    private BaseInformationRecords.CountInfo simplify(int outputComplexityLevel, BaseInformationRecords.CountInfo count) {
        if (outputComplexityLevel == 2) {
            if ((count.getGenotypeCountReverseStrand() + count.getGenotypeCountForwardStrand()) == 0) return null;
            if (args().indels && !count.getIsIndel()) {
                return null;
            }
            BaseInformationRecords.CountInfo.Builder builder = count.toBuilder();
            if (outputComplexityLevel == 2) {
                // level 2 removes all frequencies:
                builder.clearDistanceToEndOfRead();
                builder.clearDistanceToStartOfRead();
                builder.clearDistancesToReadVariationsForwardStrand();
                builder.clearDistancesToReadVariationsReverseStrand();
                builder.clearInsertSizes();
                builder.clearPairFlags();
                builder.clearTargetAlignedLengths();
                builder.clearQueryPositions();
                builder.clearQueryAlignedLengths();
                builder.clearReadIndicesForwardStrand();
                builder.clearReadIndicesReverseStrand();
                builder.clearReadMappingQualityForwardStrand();
                builder.clearReadMappingQualityReverseStrand();
                builder.clearQualityScoresForwardStrand();
                builder.clearQualityScoresReverseStrand();
                builder.clearNumVariationsInReads();
            }
            return builder.build();
        }else {
            return count;
        }
    }


    @Override
    public PrintArguments createArguments() {
        return new PrintArguments();
    }


}
