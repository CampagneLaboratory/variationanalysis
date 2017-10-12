package org.campagnelab.dl.genotype.tools;

import org.apache.commons.lang3.tuple.Pair;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.BasenameUtils;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.campagnelab.goby.util.FileExtensionHelper;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Tool to convert from SBI to SSI format.
 * 
 * @author manuele
 */
public class SBIToSSIConverter extends AbstractTool<SBIToSSIConverterArguments> {

    SequenceSegmentInformationWriter writer = null;

    SegmentList segmentList ;

    public static void main(String[] args) {
        SBIToSSIConverter tool = new SBIToSSIConverter();
        tool.parseArguments(args, "SBIToSSIConverter", tool.createArguments());
        tool.execute();
    }
    
    @Override
    public SBIToSSIConverterArguments createArguments() {
        return new SBIToSSIConverterArguments();
    }

    @Override
    public void execute() {
        if (args().inputFile.isEmpty()) {
            System.err.println("You must provide input SBI files.");
        }
        int gap = args().gap;
        try {
            if (args().ssiPrefix != null)
                writer = new SequenceSegmentInformationWriter(args().ssiPrefix);
            else
                writer = new SequenceSegmentInformationWriter(BasenameUtils.getBasename(args().inputFile,
                        FileExtensionHelper.COMPACT_SEQUENCE_BASE_INFORMATION));
            RecordReader sbiReader = new RecordReader(new File(args().inputFile).getAbsolutePath());
            BaseInformationRecords.BaseInformation sbiRecord = sbiReader.nextRecord();
            while (sbiRecord!=null) {
                manageRecord(sbiRecord, gap);
                sbiRecord = sbiReader.nextRecord();
            }
            closeOutput();
        } catch (IOException e) {
            System.err.println("Failed to parse " + args().inputFile);
            e.printStackTrace();
        }
        
    }

    private void manageRecord(BaseInformationRecords.BaseInformation record, int gap) {
        if (segmentList == null) {
            segmentList = new SegmentList(record, this.writer, null);
        } else {
            if (!this.isSameSegment(record,gap)) {
                segmentList.newSegment(record);
            } else {
                segmentList.add(record);
            }
        }
    }

    private boolean isSameSegment(BaseInformationRecords.BaseInformation record, int gap) {
         return ((record.getPosition() - segmentList.getCurrentLocation() <= gap) &&
                 (record.getReferenceId().equalsIgnoreCase(segmentList.getCurrentLastReferenceId())) );
    }

    private void closeOutput() {
        segmentList.close();
        try {
            writer.close();
        } catch (IOException e) {
            System.err.println("Failed to close the SSI file");
            e.printStackTrace();
        } finally {
            writer = null;
        }
    }



}
