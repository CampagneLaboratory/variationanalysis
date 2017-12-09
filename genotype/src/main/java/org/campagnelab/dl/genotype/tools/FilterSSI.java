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
        SBIToSSIConverter tool = new SBIToSSIConverter();
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
        try (SegmentReader ssiReader = new SegmentReader(new File(args().inputFile).getAbsolutePath());
             SegmentWriter ssiwriter = new SegmentWriter(new File(args().outputFile).getAbsolutePath());) {
            ssiReader.forEach(segmentInformation -> {
                if (this.accept(segmentInformation)) {
                    ssiwriter.writeRecord(segmentInformation);
                }
            });
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private boolean accept(SegmentInformationRecords.SegmentInformation segment) {
        if (!Objects.equals(segment.getStartPosition().getReferenceId(), this.referenceId))
            return false;
        return (segment.getStartPosition().getLocation() >= this.startLocation
                || segment.getEndPosition().getLocation() <= this.endLocation);
    }
}
