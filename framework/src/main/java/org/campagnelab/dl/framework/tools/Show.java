package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import it.unimi.dsi.io.FastBufferedReader;
import it.unimi.dsi.io.LineIterator;
import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.framework.domains.DomainDescriptor;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Iterator;
import java.util.function.Function;

/**
 * Show records of a dataset.
 */
public abstract class Show<RecordType, Args extends ShowArguments> extends AbstractTool<Args> {


    protected abstract DomainDescriptor<RecordType> domainDescriptor();

    @Override
    public void execute() {

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

        Function<RecordType, String> converter = getConverter(args().reportType.toUpperCase());
        Iterable input = domainDescriptor().getRecordIterable().apply(args().datasetFilename);
        Iterator<RecordType> reader = input.iterator();
        while (reader.hasNext()) {
            selectedIndex = getNextIndex(predictionLine, selectedIndex, index);
            RecordType next = reader.next();
            if (selectedIndex == index) {
                System.out.println(Integer.toString(index) + "\t" + converter.apply(next));
            }
            index += 1;
            if (index > args().showN) {
                break;
            }
        }
    }

    protected abstract Function<RecordType, String> getConverter(String reportType);

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
            if (line.startsWith("index")) continue;
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


}

