package org.campagnelab.dl.somatic.tools;

import org.campagnelab.dl.framework.tools.ShowArguments;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.function.Function;

/**
 * Arguments for ShowSomatic.
 */
public class SomaticShowArguments extends ShowArguments {
    @Override
    protected String defaultReportType() {
        return SomaticShowReportTypes.POSITIONS.toString();
    }

    public Function<BaseInformationRecords.BaseInformation, String> getConverter() {

        switch (SomaticShowReportTypes.valueOf(reportType)) {
            case PROTOBUFF:
                return showProtobuff;
            case COUNTS:
                return showFormattedCounts;
            case FREQ_COUNTS:
                return showFreqFormattedCounts;
            case POSITIONS:
            default:
                return showPositions;
        }
    }

    public enum SomaticShowReportTypes {
        PROTOBUFF,
        POSITIONS,
        COUNTS,
        FREQ_COUNTS
    }

    Function<BaseInformationRecords.BaseInformation, String> showFormattedCounts = new Function<BaseInformationRecords.BaseInformation, String>() {
        @Override
        public String apply(BaseInformationRecords.BaseInformation record) {
            StringBuffer formattedCounts=new StringBuffer();
            for (BaseInformationRecords.SampleInfo sample: record.getSamplesList()) {
                formattedCounts.append(sample.getFormattedCounts()+"\t");
            }
            return formattedCounts.toString();
        }
    };
    Function<BaseInformationRecords.BaseInformation, String> showFreqFormattedCounts = new Function<BaseInformationRecords.BaseInformation, String>() {
        @Override
        public String apply(BaseInformationRecords.BaseInformation record) {
            StringBuffer formattedCounts=new StringBuffer();
            for (BaseInformationRecords.SampleInfo sample: record.getSamplesList()) {
                formattedCounts.append(sample.getFormattedCounts()+"\t");
            }
            return formattedCounts.toString()+Float.toString(record.getFrequencyOfMutation());
        }
    };
    Function<BaseInformationRecords.BaseInformation, String> showPositions = new Function<BaseInformationRecords.BaseInformation, String>() {
        @Override
        public String apply(BaseInformationRecords.BaseInformation baseInformation) {
            String refId = baseInformation.hasReferenceId() ? baseInformation.getReferenceId() :
                    Integer.toString(baseInformation.getReferenceIndex());
            return String.format("%s\t%d", refId, baseInformation.getPosition());
        }
    };

    Function<BaseInformationRecords.BaseInformation, String> showProtobuff = new Function<BaseInformationRecords.BaseInformation, String>() {
        @Override
        public String apply(BaseInformationRecords.BaseInformation baseInformation) {
            return baseInformation.toString();
        }
    };
}
