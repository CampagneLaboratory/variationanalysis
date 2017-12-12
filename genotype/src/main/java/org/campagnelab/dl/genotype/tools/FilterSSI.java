package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.storage.SegmentReader;
import org.campagnelab.dl.genotype.storage.SegmentWriter;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.Objects;

/**
 * Filter SSI to extract a subset of the segments at specific positions.
 */
public class FilterSSI extends AbstractTool<FilterSSIArguments> {

    static private Logger LOG = LoggerFactory.getLogger(FilterSSI.class);
    int startLocation;
    int endLocation;
    String referenceId;

    public static void main(String[] args) {
        FilterSSI tool = new FilterSSI();
        tool.parseArguments(args, "FilterSSI", tool.createArguments());
        tool.execute();
    }

    @Override
    public FilterSSIArguments createArguments() {
        return new FilterSSIArguments();
    }

    @Override
    public void execute() {
        String[] start = this.args().startPosition.split(":", 2);
        this.referenceId = start[0];
        this.startLocation = Integer.valueOf(start[1]);
        String[] end = this.args().endPosition.split(":", 2);
        this.endLocation = Integer.valueOf(end[1]);
        if (!Objects.equals(this.referenceId, end[0]))
            throw new IllegalArgumentException("Chromosomes in the start and end positions do not match.");
        final int[] numOfFilteredBases = {0};
        try (SegmentReader ssiReader = new SegmentReader(new File(args().inputFile).getAbsolutePath());
             SegmentWriter ssiwriter = new SegmentWriter(new File(args().outputFile).getAbsolutePath());) {
            ssiReader.forEach(segmentInformation -> {
                if (this.accept(segmentInformation)) {
                    ssiwriter.writeRecord(segmentInformation);
                    numOfFilteredBases[0] += segmentInformation.getLength();
                }
            });
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Number of filtered bases: " + numOfFilteredBases[0]);

    }
    public FilterSSIArguments args() {
        return arguments;
    }

    private boolean accept(SegmentInformationRecords.SegmentInformation segment) {
        if (!Objects.equals(segment.getStartPosition().getReferenceId(), this.referenceId))
            return false;
        //condition 1: the range includes the segment
        //
        // ---startLocation---|--------------|---endLocation
        //                   seg-start      seg-end
        if ((segment.getStartPosition().getLocation() >= this.startLocation
                && segment.getEndPosition().getLocation() <= this.endLocation))
            return true;

        //condition 2: the range is fully inside the segment
        //
        // ----|---startLocation------endLocation-------|---
        //   seg-start                              seg-end
        if ((segment.getStartPosition().getLocation() <= this.startLocation
                && segment.getEndPosition().getLocation() >= this.endLocation))
            return true;

        //condition 3: the range is partially inside the segment
        //
        // ----|---startLocation------------|---endLocation----
        //   seg-start                    seg-end
        if ((segment.getStartPosition().getLocation() <= this.startLocation
                && segment.getEndPosition().getLocation() >= this.startLocation))
            //start location is in the middle of the segment
            return true;

        if ((segment.getStartPosition().getLocation() <= this.endLocation
                && segment.getEndPosition().getLocation() >= this.endLocation))
            //end location is in the middle of the segment
            return true;

        return false;
    }
}
