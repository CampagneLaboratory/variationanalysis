package org.campagnelab.dl.genotype.tools;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.floats.FloatList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSegmentsLSTM;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.genotype.mappers.NumDistinctAllelesLabelMapper;
import org.campagnelab.dl.genotype.segments.Segment;
import org.campagnelab.dl.genotype.segments.SegmentHelper;
import org.campagnelab.dl.genotype.segments.SegmentLabelMapper;
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
import java.util.Set;
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
    private static Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeaturesFunction;
    private static ThreadLocal<SegmentHelper> segmentHelper;

    // any genomic site that has strictly more indel supporting reads than the below threshold will be marked has candidateIndel.
    private int candidateIndelThreshold = 0;

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
        segmentHelper= new ThreadLocal<SegmentHelper>() {
            @Override
            protected SegmentHelper initialValue() {
                return new SegmentHelper(writer, processSegmentFunction, fillInFeaturesFunction,args().collectStatistics);
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

        processSegmentFunction = segment -> {
            segment.populateTrueGenotypes();
            BaseInformationRecords.BaseInformation previous = null;
            BaseInformationRecords.BaseInformation buildFrom = null;

            for (BaseInformationRecords.BaseInformation record : segment.recordList) {

                buildFrom = record;
                if (record.getTrueGenotype().length() > 3) {
                    Set<String> alleles = GenotypeHelper.getAlleles(record.getTrueGenotype());
                    for (String allele : alleles) {
                        if (allele.length() > 1) {
                            String insertionOrDeletion = getInsertedDeleted(allele);
                            int offset = 1;
                            for (char insertedDeleted : insertionOrDeletion.toCharArray()) {
                                previous = segment.recordList.insertAfter(previous, buildFrom, insertedDeleted, offset++);
                            }
                        }
                    }
                    // long genotypes need to be trimmed to the first base after this point. We don't do it here because
                    // it would modify the data structure currently traversed.
                }

                previous = record;
            }

            segment.recordList.forEach(record -> {
                int longestIndelLength = 0;
                for (BaseInformationRecords.SampleInfo sample : record.getSamplesList()) {
                    for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                        if (count.getIsIndel()) {
                            longestIndelLength = Math.max(longestIndelLength, count.getFromSequence().length());
                            longestIndelLength = Math.max(longestIndelLength, count.getToSequence().length());
                        }
                    }

                }
                for (int offset = 1; offset < longestIndelLength; offset++) {
                    BaseInformationRecords.BaseInformation.Builder copy = record.toBuilder();
                    //    System.out.printf("record position: %d %n",record.getPosition());
                    copy = segment.recordList.adjustCounts(copy, offset);
                    segment.insertAfter(record, copy);
                }


            });

            return segment;
        };
        SegmentLabelMapper labelMapper = new SegmentLabelMapper(args().ploidy);


        fillInFeaturesFunction = baseInformation -> {
            FloatList features = new FloatArrayList(featureMapper.numberOfFeatures());
            SegmentInformationRecords.Base.Builder builder = SegmentInformationRecords.Base.newBuilder();
            String trueGenotype = baseInformation.getTrueGenotype();
            builder.setHasCandidateIndel(hasCandidateIndel(baseInformation));
            builder.setHasTrueIndel(
                    GenotypeHelper.isIndel(baseInformation.getReferenceBase(), baseInformation.getTrueGenotype()));
            builder.setIsVariant(
                    GenotypeHelper.isVariant(baseInformation.getTrueGenotype(), baseInformation.getReferenceBase()));

            if (args().mapFeatures) {
                featureMapper.prepareToNormalize(baseInformation, 0);
                if (trueGenotype.length() > 3) {
                    //    System.out.println("Indel:" + baseInformation.getTrueGenotype());
                }
                features.clear();
                for (int featureIndex = 0; featureIndex < featureMapper.numberOfFeatures(); featureIndex++) {
                    features.add(featureMapper.produceFeature(baseInformation, featureIndex));
                }
                builder.clearFeatures();
                builder.addAllFeatures(features);
            }
            if (args().mapLabels) {
                if (trueGenotype.length() == 1) {
                    trueGenotype = trueGenotype + "|" + trueGenotype;
                }
                if (trueGenotype.length() > 3) {
                    trueGenotype = trimTrueGenotype(trueGenotype);
                }
                float[] labels = labelMapper.map(trueGenotype.replaceAll("\\|", "/"));
                builder.clearLabels();
                for (float labelValue : labels) {
                    builder.addLabels(labelValue);
                }
            }
            return builder;
        };

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

            StreamSupport.stream(sbiReader.spliterator(), args().parallel).forEach(sbiRecord -> {
                manageRecord(sbiRecord, gap);
                synchronized (pg) {
                    pg.lightUpdate();
                }
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
        // determine if a genomic site has some reads that suggest an indel.
        for (BaseInformationRecords.SampleInfo sample : baseInformation.getSamplesList()) {
            for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                if (count.getIsIndel()) {
                    if (count.getGenotypeCountForwardStrand() > candidateIndelThreshold || count.getGenotypeCountReverseStrand() > candidateIndelThreshold) {

                        return true;
                    }
                }
            }
        }
        return false;
    }

    private String trimTrueGenotype(String trueGenotype) {
        int secondBaseIndex = trueGenotype.indexOf("|");
        final String trimmed = trueGenotype.charAt(0) + "|" + trueGenotype.charAt(secondBaseIndex - 1);
        return trimmed;
    }

    private String getInsertedDeleted(String allele) {

        if (allele.charAt(1) == '-') {
            // number of deletions in allele.
            int deletionCount = 0;
            // a certain number of bases are deleted in this allele.
            for (int i = 1; i < allele.length(); i++) {
                if (allele.charAt(i) == '-') {
                    deletionCount++;
                } else {
                    break;
                }
            }
            return allele.substring(1, deletionCount);
        } else {
            return allele.substring(1, allele.length() - 2);
        }
    }

    private void manageRecord(BaseInformationRecords.BaseInformation record, int gap) {

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
    private void closeOutput() {
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
