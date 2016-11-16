package org.campagnelab.dl.somatic.tools;

import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;

import java.io.*;
import java.util.function.Function;

/**
 * Show records of a dataset.
 * TODO Show only elements with index produced by Predict (allow to chain Predict | Show)
 */
public class Show extends AbstractTool<ShowArguments> {
    public static void main(String[] args) {

        Show show = new Show();
        show.parseArguments(args, "Show", show.createArguments());
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
        Reader predReader = null;
        LineIterator predictionLine = null;
        if (args().predictionFilter != null) {
            try {
                if ("-".equals(args().predictionFilter)) {
                    predReader = new InputStreamReader(System.in);
                } else {
                    predReader = new FileReader(args().predictionFilter);
                }
            } catch (FileNotFoundException e) {
                System.err.println("Unable to open prediction filter file. -p " + args().predictionFilter);
                System.exit(1);
            }
            predictionLine = new LineIterator(new FastBufferedReader(predReader));
        }

        int index = 0;
        int selectedIndex = -1;

        Function<BaseInformationRecords.BaseInformation, String> converter;
        switch (args().reportType) {
            case PROTOBUFF:
                converter = showProtobuff;
                break;
            case POSITIONS:
            default:
                converter = showPositions;
                break;

        }
        while (reader.hasNext()) {
            selectedIndex = getNextIndex(predictionLine, selectedIndex, index);
            BaseInformationRecords.BaseInformation next = reader.next();
            if (selectedIndex == index) {
                System.out.println(Integer.toString(index) + "\t" + converter.apply(next));
            }
            index += 1;
            if (index > args().showN) {
                break;
            }
        }
    }

    /**
     * Obtain the next selectedIndex from prediction line that is strictly larger than currentIndex.
     * If predictionLine is null (option -p was not given), return currentIndex.
     *
     * @param predictionLine
     * @param selectedIndex
     * @param currentIndex
     * @return
     */
    private int getNextIndex(LineIterator predictionLine, int selectedIndex, int currentIndex) {
        if (predictionLine == null) return currentIndex;
        if (selectedIndex >= currentIndex) {
            return selectedIndex;
        }
        while (selectedIndex < currentIndex && predictionLine.hasNext()) {
// find the next selected index in the prediction file.
            MutableString line = predictionLine.next();
            int tabIndex = line.indexOf('\t');
            if (tabIndex == -1) {
                tabIndex = line.length();
            }
            if (tabIndex > 0) {
                selectedIndex = Integer.parseInt(line.subSequence(0, tabIndex).toString());
            }
        }
        return selectedIndex;

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

