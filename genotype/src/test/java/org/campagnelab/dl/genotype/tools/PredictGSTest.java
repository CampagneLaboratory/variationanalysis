package org.campagnelab.dl.genotype.tools;

import org.campagnelab.dl.genotype.predictions.SegmentGenotypePrediction;
import org.campagnelab.dl.genotype.predictions.SegmentPrediction;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by mas2182 on 12/6/17.
 */
public class PredictGSTest {

    @Test
    public void processPredictions() throws Exception {
        SegmentInformationRecords.SegmentInformation.Builder builder = SegmentInformationRecords.SegmentInformation.newBuilder();
        SegmentInformationRecords.ReferencePosition.Builder refBuilder = SegmentInformationRecords.ReferencePosition.newBuilder();
        refBuilder.setLocation(104364294);
        refBuilder.setReferenceIndex(43);
        refBuilder.setReferenceId("chr8");
        builder.setStartPosition(refBuilder.build());
        refBuilder.setLocation(104364295);
        refBuilder.setReferenceIndex(43);
        refBuilder.setReferenceId("chr8");
        builder.setEndPosition(refBuilder.build());
        SegmentInformationRecords.Sample.Builder sampleBuilder =SegmentInformationRecords.Sample.newBuilder();
        SegmentInformationRecords.Base.Builder baseBuilder = SegmentInformationRecords.Base.newBuilder();
        baseBuilder.addTrueLabel("C");
        baseBuilder.setHasCandidateIndel(true);
        baseBuilder.setHasTrueIndel(false);
        baseBuilder.setReferenceAllele("T");
        baseBuilder.setIsVariant(true);
        baseBuilder.setFormattedCounts(" G=0  A=0  C=43  T=0  N=0  (from: T)");
        baseBuilder.setPrePostProcessingGenotype("");
        baseBuilder.setLocation(104364294);
        sampleBuilder.addBase(baseBuilder.build());
        baseBuilder = SegmentInformationRecords.Base.newBuilder();
        baseBuilder.addTrueLabel("G");
        baseBuilder.addTrueLabel("T");
        baseBuilder.setHasCandidateIndel(true);
        baseBuilder.setHasTrueIndel(false);
        baseBuilder.setReferenceAllele("G");
        baseBuilder.setIsVariant(true);
        baseBuilder.setFormattedCounts(" G=137  C=0  N=0  T=123  A=0  (from: G)");
        baseBuilder.setPrePostProcessingGenotype("");
        baseBuilder.setLocation(104364295);
        sampleBuilder.addBase(baseBuilder.build());
        builder.addSample(sampleBuilder.build());
        builder.setLength(sampleBuilder.getBaseCount());
        SegmentGenotypePrediction segmentGenotypePrediction = new SegmentGenotypePrediction();
        segmentGenotypePrediction.predictedGenotypes = new String[] {"CC", "G-"};
        segmentGenotypePrediction.length = builder.getLength();
        segmentGenotypePrediction.probabilities = new float[] {0.123F, 0.333F};
        segmentGenotypePrediction.trueGenotypes = new String[] {"",""};
        SegmentPrediction fullPred = new SegmentPrediction(builder.getStartPosition(),
                builder.getEndPosition(),segmentGenotypePrediction);

        PredictGS tool = new PredictGS();
        String[] args = new String[]{"-i", "1234-MM-test", "-m", "fromTest"};
        tool.parseArguments(args,"predictGS",tool.createArguments());
        tool.writeHeader(null);
        tool.processAggregatedPrediction(null,builder.build(),fullPred);
        tool.reportStatistics("");
    }

}