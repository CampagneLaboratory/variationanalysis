package org.campagnelab.dl.genotype.mappers;

import com.google.protobuf.TextFormat;
import org.apache.commons.lang.ArrayUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by rct66 on 12/5/16.
 */
public class GenotypeLabelMapperTest {
    @Test
    public void mapFeatures() throws Exception {
        int index=0;

        for (String record : records) {
            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);


            for (int i = 0; i < 10; i++){
                GenotypeLabelsMapper calculator = new GenotypeLabelsMapper(i);
                INDArray inputs = Nd4j.zeros(ArrayUtils.addAll(new int[]{1}, calculator.dimensions().dimensions));
                calculator.mapLabels(builder.build(),inputs,0);
                assertEquals(expectedLabels[index][i], inputs.toString());


            }
            index++;

        }
    }

    String[] records = {
            "reference_index: 21\n" +
                    "position: 9413944\n" +
                    "mutated: true\n" +
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
                    "    genotypeCountForwardStrand: 4\n" +
                    "    genotypeCountReverseStrand: 5\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 127\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 138\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 133\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 131\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 106\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 115\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 121\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: false\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 12\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 33\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 17\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 3\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 7\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 30\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -848\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -578\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -689\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 331\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -573\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 363\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 489\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 428\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -652\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 34\n" +
                    "    genotypeCountReverseStrand: 33\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 14\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 9\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 16\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 250\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 247\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 246\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 238\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 234\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 232\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 231\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 223\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 217\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 192\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 191\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 164\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 160\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 155\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 146\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 145\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 113\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 100\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 93\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 92\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 86\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 83\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 61\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 59\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 46\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 43\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 34\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 17\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 5\n" +
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
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 25\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 26\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 30\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 31\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 46\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 63\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 65\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 88\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 103\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 108\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 135\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 142\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 148\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 149\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 159\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 175\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 192\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 194\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 196\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 198\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 200\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 201\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 205\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 206\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 223\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 228\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 236\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 239\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 30\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 30\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 3\n" +
                    "      frequency: 28\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 6\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 8\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 5\n" +
                    "      frequency: 3\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 4\n" +
                    "      frequency: 9\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 7\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 2\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 453\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 566\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 643\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 386\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -419\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 665\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -430\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 896\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 430\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 364\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -469\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -435\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 526\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -446\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -537\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -474\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 838\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 578\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 711\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -644\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -568\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -732\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -580\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 441\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 583\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 345\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -486\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 412\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 366\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -361\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -364\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 467\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 463\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -386\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -600\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -465\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 450\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 357\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 501\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 421\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -470\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 397\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -440\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 521\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -345\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 552\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -560\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -331\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -640\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -515\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -453\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 413\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -467\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 455\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 407\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 488\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -366\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -540\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 405\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 398\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -536\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -615\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: false\n" +
                    "    fromSequence: \"T\"\n" +
                    "    toSequence: \"G\"\n" +
                    "    genotypeCountForwardStrand: 4\n" +
                    "    genotypeCountReverseStrand: 1\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 201\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 38\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 14\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 8\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 34\n" +
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
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 11\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 0\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 19\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 37\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 8\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: -754\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 518\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 737\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 704\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 760\n" +
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
                    "  formattedCounts: \"sample: 0 counts A=0 T=9 C=67 G=5 N=0 FB=0 indels={ [] }\\n\"\n" +
                    "}"
    };
    String[][] expectedLabels = {{
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[1.00, 0.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]",
            "[0.00, 1.00]"}
            };
}