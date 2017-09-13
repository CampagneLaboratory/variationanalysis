package org.campagnelab.dl.genotype.mappers;

import com.google.protobuf.TextFormat;
import org.apache.commons.lang.ArrayUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import static org.junit.Assert.assertEquals;

/**
 * Created by rct66 on 12/5/16.
 */
public class CombinedGenotypeLabelMapperTest {
    @Test
    public void mapFeatures() throws Exception {
        int index=0;

        for (String record : records) {
            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);


            for (int i = 0; i < 1; i++){
                CombinedLabelsMapper calculator = new CombinedLabelsMapper();
                INDArray inputs = Nd4j.zeros(ArrayUtils.addAll(new int[]{1}, calculator.dimensions().dimensions));
                calculator.prepareToNormalize(builder.build(),0);
                calculator.mapLabels(builder.build(),inputs,0);
                assertEquals(expectedLabels[index][i], inputs.toString());
            }
            index++;

        }
    }

    String[] records = {
            "   reference_index: 21\n" +
                    "position: 29382312\n" +
                    "mutated: false\n" +
                    "referenceBase: \"T\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    isCalled: false\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 14\n" +
                    "    genotypeCountReverseStrand: 24\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 14\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 24\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 1\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 4\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 35\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 109\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 115\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 162\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 163\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 165\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 177\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 180\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 198\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 226\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 246\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 12\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 16\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 44\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 56\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 74\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 80\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 90\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 91\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 99\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 109\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 114\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 126\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 142\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 166\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 169\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 191\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 210\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 216\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 238\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 248\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 15\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 0\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 1\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 3\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 4\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 6\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 10\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 12\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 14\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 21\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 24\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 32\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 35\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -781\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -746\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -635\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -631\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -614\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -603\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -562\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -542\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -533\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -532\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -531\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -528\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -510\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -491\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -489\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -463\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -462\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -456\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -427\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -426\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -418\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -399\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -395\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -303\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 411\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 427\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 447\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 477\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 486\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 525\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 551\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 556\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 564\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 579\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 594\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 597\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 678\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 719\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    isCalled: false\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"G\"\n" +
                    "    genotypeCountForwardStrand: 6\n" +
                    "    genotypeCountReverseStrand: 6\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 5\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 24\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 91\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 92\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 179\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 241\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 14\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 44\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 132\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 200\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 202\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 212\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 3\n" +
                    "      frequency: 7\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 14\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 21\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 36\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -755\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -668\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -587\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -550\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -471\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -470\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 486\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 497\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 505\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 533\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 710\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 736\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    isCalled: false\n" +
                    "  }\n" +
                    "  isTumor: true\n" +
                    "  formattedCounts: \"sample: 0 counts A=0 T=38 C=0 G=12 N=0 FB=0 indels={ null }\\n\"\n" +
                    "  isVariant: true\n" +
                    "}\n" +
                    "trueGenotype: \"T/G\"\n" +
                    "reference_id: \"chr21\"\n" +
                    "genomicSequenceContext: \"TTGCATGTTTTTTACTTTTTT\"\n"
            ,
            "reference_index: 21\n" +
                    "position: 9872928\n" +
                    "mutated: false\n" +
                    "referenceBase: \"G\"\n" +
                    "samples {\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"G\"\n" +
                    "    toSequence: \"A\"\n" +
                    "    genotypeCountForwardStrand: 8\n" +
                    "    genotypeCountReverseStrand: 8\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 35\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 73\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 80\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 108\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 161\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 236\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 240\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 31\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 48\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 107\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 117\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 223\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 235\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 243\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 1\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 8\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 25\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -664\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -561\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -555\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -493\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -469\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -397\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -344\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -301\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 301\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 408\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 444\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 469\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 511\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 544\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 568\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 627\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"G\"\n" +
                    "    toSequence: \"T\"\n" +
                    "    genotypeCountForwardStrand: 1\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 178\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: false\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 546\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"G\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 2\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 133\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 189\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: false\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 53\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 489\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 621\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"G\"\n" +
                    "    toSequence: \"G\"\n" +
                    "    genotypeCountForwardStrand: 93\n" +
                    "    genotypeCountReverseStrand: 103\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 93\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 103\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 17\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 28\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 30\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 36\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 37\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 41\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 43\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 44\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 51\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 62\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 63\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 64\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 73\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 76\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 81\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 82\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 87\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 89\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 90\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 91\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 95\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 98\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 101\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 107\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 108\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 109\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 115\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 117\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 118\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 123\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 124\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 127\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 128\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 129\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 131\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 132\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 136\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 140\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 147\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 152\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 155\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 158\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 161\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 163\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 166\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 170\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 171\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 176\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 177\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 181\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 189\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 191\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 193\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 201\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 202\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 203\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 208\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 211\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 213\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 220\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 222\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 223\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 229\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 233\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 234\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 235\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 239\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 242\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 244\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 245\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 246\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 250\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 6\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 10\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 11\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 18\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 21\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 24\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 25\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 42\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 44\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 45\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 48\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 50\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 55\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 66\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 69\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 72\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 74\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 77\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 82\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 86\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 88\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 89\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 92\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 95\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 99\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 103\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 104\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 105\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 107\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 109\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 110\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 113\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 115\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 120\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 121\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 122\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 123\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 128\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 130\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 133\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 138\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 140\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 150\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 152\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 153\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 156\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 157\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 163\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 164\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 166\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 169\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 171\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 172\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 175\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 176\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 191\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 192\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 196\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 197\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 199\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 200\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 205\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 211\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 214\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 221\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 225\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 227\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 228\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 230\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 232\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 234\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 236\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 240\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 245\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 246\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 247\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 248\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 251\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 7\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 17\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 20\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 15\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 43\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 46\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 67\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 17\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 20\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 37\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 76\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 0\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 1\n" +
                    "      frequency: 51\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 33\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 3\n" +
                    "      frequency: 22\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 4\n" +
                    "      frequency: 23\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 6\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 7\n" +
                    "      frequency: 7\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 8\n" +
                    "      frequency: 15\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 9\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 10\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 11\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 12\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 14\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 19\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 20\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 31\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 37\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 42\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -838\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -748\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -745\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -686\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -675\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -673\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -670\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -668\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -664\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -658\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -639\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -613\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -604\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -602\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -597\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -595\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -591\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -587\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -583\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -575\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -571\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -570\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -555\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -553\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -551\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -546\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -542\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -534\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -524\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -519\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -507\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -506\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -504\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -496\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -495\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -493\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -491\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -485\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -483\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -477\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -475\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -474\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -468\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -466\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -465\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -459\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -456\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -455\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -449\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -446\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -442\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -441\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -440\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -437\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -435\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -433\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -432\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -431\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -430\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -429\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -428\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -427\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -424\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -421\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -417\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -414\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -409\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -408\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -402\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -400\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -398\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -396\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -395\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -394\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -392\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -391\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -390\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -383\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -382\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -380\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -376\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -374\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -372\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -364\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -359\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -334\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -313\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -311\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -305\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -299\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -287\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -242\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 242\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 287\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 311\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 313\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 345\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 359\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 364\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 367\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 372\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 374\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 382\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 386\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 390\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 395\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 400\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 401\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 402\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 403\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 406\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 407\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 417\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 419\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 421\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 429\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 430\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 432\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 436\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 437\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 440\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 442\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 443\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 446\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 447\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 456\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 457\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 464\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 465\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 468\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 469\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 471\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 474\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 478\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 480\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 485\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 486\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 487\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 489\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 492\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 493\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 495\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 504\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 516\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 526\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 528\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 548\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 559\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 561\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 564\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 565\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 574\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 578\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 583\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 585\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 587\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 591\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 592\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 593\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 599\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 606\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 609\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 619\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 635\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 636\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 637\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 652\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 670\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 672\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 674\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 694\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 703\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 716\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 717\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 718\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 735\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 749\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 807\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"G\"\n" +
                    "    toSequence: \"N\"\n" +
                    "    genotypeCountForwardStrand: 0\n" +
                    "    genotypeCountReverseStrand: 0\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 0\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 154\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: false\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 531\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  isTumor: true\n" +
                    "  formattedCounts: \"sample: 0 counts A=16 T=1 C=2 G=196 N=0 FB=0 indels={ null }\\n\"\n" +
                    "  isVariant: true\n" +
                    "}\n" +
                    "trueGenotype: \"G/A\"\n" +
                    "reference_id: \"chr21\"\n" +
                    "genomicSequenceContext: \"TCAACATTAGGCTCAGGTGGA\"\n" +
                    "\n"


    };
    String[][] expectedLabels = {
            {"[0.00,  1.00,  0.00,  0.00]"},
            {"[0.00,  1.00,  0.00,  0.00]"}

    };
}