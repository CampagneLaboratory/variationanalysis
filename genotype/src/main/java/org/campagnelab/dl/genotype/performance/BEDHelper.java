package org.campagnelab.dl.genotype.performance;

import org.apache.commons.io.IOUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Helper to write bed file for observed regions, as well as TP, TN, FP, FN.
 * Created by fac2003 on 2/20/17.
 */
public class BEDHelper {
    private PrintWriter positionWriter;
    private PrintWriter fpWriter;
    private PrintWriter fnWriter;
    private PrintWriter tpWriter;
    private PrintWriter tnWriter;

    public BEDHelper(String basename) throws IOException {
        positionWriter = new PrintWriter(new FileWriter(basename + "-observed-regions.bed"));
        fpWriter = new PrintWriter(new FileWriter(basename + "-fp.bed"));
        fnWriter = new PrintWriter(new FileWriter(basename + "-fn.bed"));
        tpWriter = new PrintWriter(new FileWriter(basename + "-tp.bed"));
        tnWriter = new PrintWriter(new FileWriter(basename + "-tn.bed"));
    }

    public void close() {
        IOUtils.closeQuietly(positionWriter);
        IOUtils.closeQuietly(fpWriter);
        IOUtils.closeQuietly(fnWriter);
        IOUtils.closeQuietly(tpWriter);
        IOUtils.closeQuietly(tnWriter);
    }

    public void add(String referenceId, int start, int end, int index, StatsAccumulator stats) {
        positionWriter.printf("%s\t%d\t%d\t%d\n", referenceId, start, end, index);
        if (stats.observedWasFP()) {
            fpWriter.printf("%s\t%d\t%d\t%d\n", referenceId, start, end, index);
        }
        if (stats.observedWasFN()) {
            fnWriter.printf("%s\t%d\t%d\t%d\n", referenceId, start, end, index);
        }
        if (stats.observedWasTP()) {
            tpWriter.printf("%s\t%d\t%d\t%d\n", referenceId, start, end, index);
        }
        if (stats.observedWasTN()) {
            tnWriter.printf("%s\t%d\t%d\t%d\n", referenceId, start, end, index);
        }


    }
}
