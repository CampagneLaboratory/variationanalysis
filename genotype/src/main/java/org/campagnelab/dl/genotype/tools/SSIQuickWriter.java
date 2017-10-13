package org.campagnelab.dl.genotype.tools;

import com.google.common.primitives.Floats;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.storage.SegmentWriter;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Tool to quickly write out an SSI file for a simple learning task.
 * Write out a protobuf SSI file with --num-segments SegmentInformations, with one Sample for each
 * SegmentInformation, --bases-per-segment Bases for each Sample, and --features-per-base features for each Base,
 * which are randomly generated integers between 0 and --upper-bound.
 * For the set of features for each base, the labels correspond to
 * - L[0] = 1 if num_features_even is even, 0 otherwise
 * - L[1] = 1 if num_features_odd is odd, 0 otherwise
 * - L[2] = 1 if num_features_greater_than_average_of_features >= (features.length / 2)
 *
 * @author joshuacohen
 */

public class SSIQuickWriter extends AbstractTool<SSIQuickWriterArguments> {
    public static void main(String[] args) {
        SSIQuickWriter tool = new SSIQuickWriter();
        tool.parseArguments(args, "SSIQuickWriter", tool.createArguments());
        tool.execute();
    }

    @Override
    public SSIQuickWriterArguments createArguments() {
        return new SSIQuickWriterArguments();
    }

    @Override
    public void execute() {
        try {
            SegmentWriter segmentWriter = new SegmentWriter(args().outputFile);
            Random rng = new XoRoShiRo128PlusRandom();
            System.out.println("Num segments to write: " + args().numSegments);
            System.out.println("Bases per sequence: " + args().basesPerSegment);
            System.out.println("Features per base: " + args().featuresPerBase);
            System.out.println("Upper bound: " + args().upperBound);
            for (int i = 0; i < args().numSegments; i++) {
                SegmentInformationRecords.Sample.Builder sampleBuilder = SegmentInformationRecords.Sample.newBuilder();
                for (int j = 0; j < args().basesPerSegment; j++) {
                    float[] features = new float[args().featuresPerBase];
                    int totalEven = 0;
                    for (int k = 0; k < args().featuresPerBase; k++) {
                        features[k] = rng.nextInt(args().upperBound);
                        if (features[k] % 2 == 0) {
                            totalEven++;
                        }
                    }
                    int totalOdd = features.length - totalEven;
                    float[] labels = new float[3];
                    labels[0] = totalEven % 2 == 0 ? 1 : 0;
                    labels[1] = totalOdd % 2 == 1 ? 1 : 0;
                    double average = IntStream.range(0, features.length).mapToDouble(f -> features[f]).sum() / features.length;
                    long greaterThanAverage = IntStream.range(0, features.length).mapToDouble(f -> features[f]).filter(f -> f >= average).count();
                    labels[2] = greaterThanAverage >= features.length / 2 ? 1 : 0;
                    sampleBuilder.addBase(SegmentInformationRecords.Base.newBuilder()
                            .addAllFeatures(Floats.asList(features))
                            .addAllLabels(Floats.asList(labels))
                            .build()
                    );
                }
                segmentWriter.writeRecord(SegmentInformationRecords.SegmentInformation.newBuilder()
                        .addSample(sampleBuilder.build())
                        .build()
                );
            }
            segmentWriter.close();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }
}
