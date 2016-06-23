package org.campagnelab.dl.varanalysis.learning.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.learning.mappers.SimpleFeatureCalculator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.*;

/**
 * Created by fac2003 on 5/27/16.
 */
public class SimpleFeatureCalculatorTest {
    @Test
    public void mapFeatures() throws Exception {
int index=0;

        for (String record : records) {
           SimpleFeatureCalculator calculator=new SimpleFeatureCalculator();

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
            "reference_index: 18\n" +
            "position: 17214616\n" +
            "mutated: true\n" +
            "mutatedBase: \"C\"\n"+
            "indexOfMutatedBase: 2\n"+
            "referenceBase: \"A\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "}\n" +
            "samples {\n" +
            "  isTumor: true\n"+
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 50\n" +
            "    genotypeCountReverseStrand: 50\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 25\n" +
            "    genotypeCountReverseStrand: 25\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"A\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "  }\n" +
            "} " ,
            "reference_index: 18\n" +
                    "position: 17214616\n" +
                    "mutated: true\n" +
                    "mutatedBase: \"C\"\n"+
                    "indexOfMutatedBase: 2\n"+
                    "referenceBase: \"A\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 5\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 4\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 3\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"G\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 2\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 1\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "}\n" +
                    "samples {\n" +
                    "  isTumor: true\n"+
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 1\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 2\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 3\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"G\"\n" +
                    "    genotypeCountForwardStrand: 4\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"A\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 5\n" +
                    "  }\n" +
                    "} "
                    ,
                    "reference_index: 0\n" +
                            "position: 20913\n" +
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
            "[0.12, 0.12, 0.12, 0.12, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12, 0.12, 0.06, 0.06, 0.00, 0.00, 0.06, 0.06, 0.00, 0.00]",
            "[0.00, 0.17, 0.13, 0.00, 0.00, 0.10, 0.07, 0.00, 0.00, 0.03, 0.03, 0.00, 0.00, 0.07, 0.10, 0.00, 0.00, 0.13, 0.17, 0.00]",
            "[0.26, 0.20, 0.12, 0.12, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.14, 0.02, 0.02, 0.00, 0.03, 0.00, 0.00, 0.00, 0.00]"};
}