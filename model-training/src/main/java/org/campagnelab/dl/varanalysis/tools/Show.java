package org.campagnelab.dl.varanalysis.tools;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;

import java.io.IOException;
import java.util.function.Function;

/**
 * Show records of a dataset.
 * TODO Show only elements with index produced by Predict (allow to chain Predict | Show)
 */
public class Show extends AbstractTool<ShowArguments> {
    public static void main(String[] args) {

        Show show = new Show();
        show.parseArguments(args, "Predict", show.createArguments());
        show.execute();
    }

    @Override
    public ShowArguments createArguments() {
        return new ShowArguments();
    }


    @Override
    public void execute() {
        SequenceBaseInformationReader reader = null;

        try {
            reader = new SequenceBaseInformationReader(args().datasetFilename);
        } catch (IOException e) {
            System.err.println("Unable to load input dataset: " + args().datasetFilename);
            e.printStackTrace();
        }
        if (reader == null) {
            System.err.println("Unable to create reader for input dataset.");
        }
        int index = 0;
        Function<BaseInformationRecords.BaseInformation, String> converter=showPositions;
        while (reader.hasNext()) {
            BaseInformationRecords.BaseInformation next = reader.next();
            System.out.println(converter.apply(next));
            index += 1;
            if (index > args().showN) {
                break;
            }
        }
    }

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
