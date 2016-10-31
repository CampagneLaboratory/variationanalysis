package org.campagnelab.dl.varanalysis.learning.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FeatureMapperV9;
import org.campagnelab.dl.model.utils.mappers.OneHotBaseMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.assertEquals;

/**
 * Created by fac2003 on 5/27/16.
 */
public class OneHotBaseMapperTest {
    @Test
    public void mapFeatures() throws Exception {
        int index=0;

        for (String record : records) {
           FeatureMapper calculator=new OneHotBaseMapper(0, BaseInformationRecords.BaseInformationOrBuilder::getGenomicSequenceContext);

            INDArray inputs = Nd4j.zeros(1, calculator.numberOfFeatures());

            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);

            calculator.prepareToNormalize(builder.build(),0);
            calculator.mapFeatures(builder.build(),inputs,0);


            assertEquals(expectedFeatures[index], inputs.toString());
            index++;
        }
    }

    String[] records = {
                    "reference_index: 0\n" +
                            "position: 20913\n" +
                            "genomicSequenceContext: \"GTTTTTTTTTTTTTTTTTTTT\"\n" +
                            "mutated: false\n" +
                            "referenceBase: \"A\"\n" +
                            "samples {\n" +
                            "  isTumor: false\n"+
                            "  counts {\n" +
                            "    matchesReference: true\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"A\"\n" +
                            "    genotypeCountForwardStrand: 13\n" +
                            "    genotypeCountReverseStrand: 17\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 7\n" +
                            "    }\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 6\n" +
                            "    }\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 6\n" +
                            "    }\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 11\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 84\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 83\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 72\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 71\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 59\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 56\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 48\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 42\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 41\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 30\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 25\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 10\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 42\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 51\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 52\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 58\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 61\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 62\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 63\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 64\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 65\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 66\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 80\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 96\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 101\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"T\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"C\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 3\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 76\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 85\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 91\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"G\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"N\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 7\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 55\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 28\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 40\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 45\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 86\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 91\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"A---GGAG\"\n" +
                            "    genotypeCountForwardStrand: 8\n" +
                            "    genotypeCountReverseStrand: 8\n" +
                            "    isIndel: true\n" +
                            "  }\n" +
                            "  formattedCounts: \"sample: 0 counts A=30 T=0 C=3 G=0 N=0 FB=0 indels={ [indel count=8 A GGAGGAG/---GGAG  20913-20921 filtered=false] }\\n\"\n" +
                            "}\n" +
                            "samples {\n" +
                            "  isTumor: true\n"+
                            "  counts {\n" +
                            "    matchesReference: true\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"A\"\n" +
                            "    genotypeCountForwardStrand: 9\n" +
                            "    genotypeCountReverseStrand: 3\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 5\n" +
                            "    }\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 4\n" +
                            "    }\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 3\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 94\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 70\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 65\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 46\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 34\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 21\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 8\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 4\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 1\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 33\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 96\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"T\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"C\"\n" +
                            "    genotypeCountForwardStrand: 2\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresForwardStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "    readIndicesForwardStrand {\n" +
                            "      number: 47\n" +
                            "      frequency: 2\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"G\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"N\"\n" +
                            "    genotypeCountForwardStrand: 0\n" +
                            "    genotypeCountReverseStrand: 0\n" +
                            "    isIndel: false\n" +
                            "    qualityScoresReverseStrand {\n" +
                            "      number: 0\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "    readIndicesReverseStrand {\n" +
                            "      number: 90\n" +
                            "      frequency: 1\n" +
                            "    }\n" +
                            "  }\n" +
                            "  counts {\n" +
                            "    matchesReference: false\n" +
                            "    fromSequence: \"AGGAGGAG\"\n" +
                            "    toSequence: \"A---GGAG\"\n" +
                            "    genotypeCountForwardStrand: 1\n" +
                            "    genotypeCountReverseStrand: 1\n" +
                            "    isIndel: true\n" +
                            "  }\n" +
                            "  formattedCounts: \"sample: 1 counts A=12 T=0 C=2 G=0 N=0 FB=0 indels={ [indel count=1 A GGAGGAG/---GGAG  20913-20921 filtered=false] }\\n\"\n" +
                            "}"
    };
    String[] expectedFeatures = {
            "[0.00, 0.00, 0.00, 1.00, 0.00, 0.00]"};
}