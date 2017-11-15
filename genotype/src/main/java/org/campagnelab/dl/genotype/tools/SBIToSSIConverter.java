package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.genotype.mappers.NumDistinctAllelesLabelMapper;
import org.campagnelab.dl.genotype.segments.*;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.campagnelab.goby.baseinfo.BasenameUtils;
import org.campagnelab.goby.baseinfo.SequenceSegmentInformationWriter;
import org.campagnelab.goby.util.FileExtensionHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.Properties;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.StreamSupport;


/**
 * Tool to convert from SBI to SSI format.
 *
 * @author manuele
 */
public class SBIToSSIConverter extends AbstractTool<SBIToSSIConverterArguments> {

    static private Logger LOG = LoggerFactory.getLogger(SBIToSSIConverter.class);
    private static SequenceSegmentInformationWriter writer;
    private static Function<Segment, Segment> processSegmentFunction;
    private static FillInFeaturesFunction fillInFeaturesFunction;
    private static ThreadLocal<SegmentHelper> segmentHelper;

    // any genomic site that has strictly more indel supporting reads than the below threshold will be marked has candidateIndel.
    private static int candidateIndelThreshold = 0;

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
        Consumer<SegmentInformationRecords.SegmentInformation> segmentConsumer = segment -> {

            writer.appendEntry(segment);

        };

        segmentHelper = new ThreadLocal<SegmentHelper>() {
            @Override
            protected SegmentHelper initialValue() {
                SegmentHelper helper = new SegmentHelper(processSegmentFunction, fillInFeaturesFunction, segmentConsumer,
                        args().getStrategy(),
                        args().collectStatistics);
                helper.setSamplingRate(args().samplingRate);
                return helper;
            }
        };
        int gap = args().gap;
        GenotypeDomainDescriptor domainDescriptor = null;
        try {
            RecordReader sbiReader = new RecordReader(new File(args().inputFile).getAbsolutePath());

            Properties sbiProperties = sbiReader.getProperties();
            Properties domainProperties = new Properties();
            domainProperties.put("net.architecture.classname", GenotypeSegmentsLSTM.class.getCanonicalName());
            domainProperties.put(NumDistinctAllelesLabelMapper.PLOIDY_PROPERTY, Integer.toString(args().ploidy));
            domainProperties.put("input.featureMapper", args().featureMapperClassName);
            domainProperties.put("genomicContextLength", "1");
            domainProperties.put("indelSequenceLength", "1");

            domainDescriptor = new GenotypeDomainDescriptor(domainProperties, sbiProperties);
        } catch (IOException e) {
            System.err.println("Unable to initialized genotype domain descriptor.");
            e.printStackTrace();
        }


        FeatureMapper featureMapper = domainDescriptor.getFeatureMapper("input");
        //LabelMapper labelMapper=domainDescriptor.getFeatureMapper("input");
        if (args().snpOnly) {
            processSegmentFunction = new SnpOnlyPostProcessSegmentFunction();
        } else {
            processSegmentFunction = new WithIndelsPostProcessSegmentFunction();

        }
        SegmentLabelMapper labelMapper = new SegmentLabelMapper(args().ploidy);


        fillInFeaturesFunction = new MyFillInFeaturesFunction(featureMapper, labelMapper, arguments);

        try {
            if (args().ssiPrefix != null)
                writer = new SequenceSegmentInformationWriter(args().ssiPrefix);
            else
                writer = new SequenceSegmentInformationWriter(BasenameUtils.getBasename(args().inputFile,
                        FileExtensionHelper.COMPACT_SEQUENCE_BASE_INFORMATION));
            Properties props = new Properties();
            labelMapper.writeMap(props);
            writer.appendProperties(props);
            RecordReader sbiReader = new RecordReader(new File(args().inputFile).getAbsolutePath());
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.displayFreeMemory = true;
            pg.expectedUpdates = sbiReader.getTotalRecords();
            pg.itemsName = "records";
            pg.start();
            final int[] totalRecords = {0};

            StreamSupport.stream(sbiReader.spliterator(), args().parallel).limit(args().readN).forEach(sbiRecord -> {
                try {
                    manageRecord(sbiRecord, gap);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                pg.lightUpdate();
                totalRecords[0]++;
            });

            System.out.printf("Total record managed: %d %n", totalRecords[0]);
            pg.stop();
            closeOutput();
        } catch (IOException e) {
            System.err.println("Failed to parse " + args().inputFile);
            e.printStackTrace();
        }

    }

    private boolean hasCandidateIndel(BaseInformationRecords.BaseInformation baseInformation) {
        return SegmentUtil.hasCandidateIndel(baseInformation, candidateIndelThreshold);
    }

    private String trimTrueGenotype(String trueGenotype) {
        int secondBaseIndex = trueGenotype.indexOf("|");
        final String trimmed = trueGenotype.charAt(0) + "|" + trueGenotype.charAt(secondBaseIndex - 1);
        return trimmed;
    }

    private void manageRecord(BaseInformationRecords.BaseInformation record, int gap) throws IOException {

        if (this.isValid(record)) {
            final SegmentHelper segmentHelper = SBIToSSIConverter.segmentHelper.get();
            if (!this.isSameSegment(record, gap)) {
                segmentHelper.newSegment(record);
            } else {
                segmentHelper.add(record);
            }
        }
    }

    /**
     * Checks if this record should be considered.
     *
     * @param record
     * @return
     */
    private boolean isValid(BaseInformationRecords.BaseInformation record) {
        final int[] sum = {0};
        record.getSamplesList().forEach(sample -> {
                    sample.getCountsList().forEach(counts -> {
                        sum[0] += counts.getGenotypeCountForwardStrand();
                        sum[0] += counts.getGenotypeCountReverseStrand();
                    });
                }
        );
        //System.out.println("Sum for the counts is " + sum[0]);
        return (sum[0] > 0);
    }

    /**
     * Checks if the record belongs to the current segment.
     *
     * @param record
     * @param gap
     * @return
     */
    private boolean isSameSegment(BaseInformationRecords.BaseInformation record, int gap) {
        if (segmentHelper == null) {
            return false;
        }
        final boolean valid = (record.getPosition() - segmentHelper.get().getCurrentLocation() <= gap) &&
                record.getReferenceIndex() == segmentHelper.get().getCurrentReferenceIndex();
    /*if (!valid) {
        System.out.printf("not valid, actual gap=%d%n",record.getPosition() - segmentList.getCurrentLocation());
    }*/
        return valid;

    }

    /**
     * Closes the list and serializes the output SSI.
     */
    private void closeOutput() throws IOException {
        segmentHelper.get().close();
        try {
            writer.close();
        } catch (IOException e) {
            System.err.println("Failed to close the SSI file");
            e.printStackTrace();
        } finally {
            writer = null;
        }
        segmentHelper.get().printStats();
    }


}
