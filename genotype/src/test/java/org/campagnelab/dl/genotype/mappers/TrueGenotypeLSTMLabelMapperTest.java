package org.campagnelab.dl.genotype.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Properties;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 2/7/17.
 */
public class TrueGenotypeLSTMLabelMapperTest {

    @Test
    public void testRecord1() throws Exception {
        String expectedLabel = "[[[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask = "[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00]";
        TrueGenotypeLSTMLabelMapper mapper = new TrueGenotypeLSTMLabelMapper(30,0);
        MappedDimensions dim = mapper.dimensions();
        INDArray labels = Nd4j.zeros(1, dim.numElements(1), dim.numElements(2));
        INDArray mask = Nd4j.zeros(1, dim.numElements(2));
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(records[0], builder);
        final BaseInformationRecords.BaseInformation recordObject = builder.build();
        mapper.prepareToNormalize(recordObject, 0);
        mapper.mapLabels(recordObject, labels, 0);
        mapper.maskLabels(recordObject, mask, 0);
        assertEquals(expectedLabel, labels.toString());
        assertEquals(expectedMask, mask.toString());
    }

    @Test
    public void testRecord2() throws Exception {
        String expectedLabel = "[[[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00],  \n" +
                "  [1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask = "[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  0.00]";
        TrueGenotypeLSTMLabelMapper mapper = new TrueGenotypeLSTMLabelMapper(30,0);
        MappedDimensions dim = mapper.dimensions();
        INDArray labels = Nd4j.zeros(1, dim.numElements(1), dim.numElements(2));
        INDArray mask = Nd4j.zeros(1, dim.numElements(2));
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(records[1], builder);
        final BaseInformationRecords.BaseInformation recordObject = builder.build();
        mapper.prepareToNormalize(recordObject, 0);
        mapper.mapLabels(recordObject, labels, 0);
        mapper.maskLabels(recordObject, mask, 0);
        assertEquals(expectedLabel, labels.toString());
        assertEquals(expectedMask, mask.toString());
    }

    private static final String RECORD_BASE = "reference_index: 21\n" +
            "position: 45944850\n" +
            "mutated: false\n" +
            "referenceBase: \"C\"\n" +
            "samples {\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"CAGTCAGTCAGTCAGTCAGTCAGT\"\n" +
            "    toSequence: \"C----AGTCAGTCAGTCAGTCAGT\"\n" +
            "    genotypeCountForwardStrand: 19\n" +
            "    genotypeCountReverseStrand: 19\n" +
            "    isIndel: true\n" +
            "    isCalled: true\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"C\"\n" +
            "    toSequence: \"A\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "    isCalled: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"C\"\n" +
            "    toSequence: \"T\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "    isCalled: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: true\n" +
            "    fromSequence: \"C\"\n" +
            "    toSequence: \"C\"\n" +
            "    genotypeCountForwardStrand: 12\n" +
            "    genotypeCountReverseStrand: 13\n" +
            "    isIndel: false\n" +
            "    qualityScoresForwardStrand {\n" +
            "      number: 40\n" +
            "      frequency: 12\n" +
            "    }\n" +
            "    qualityScoresReverseStrand {\n" +
            "      number: 40\n" +
            "      frequency: 13\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 10\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 28\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 52\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 75\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 76\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 111\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 133\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 138\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 139\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 141\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 198\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 251\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 6\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 12\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 14\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 26\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 52\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 66\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 102\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 161\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 199\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 222\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 236\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 244\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 250\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    isCalled: true\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 29\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 50\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 60\n" +
            "      frequency: 10\n" +
            "    }\n" +
            "    readMappingQualityReverseStrand {\n" +
            "      number: 29\n" +
            "      frequency: 5\n" +
            "    }\n" +
            "    readMappingQualityReverseStrand {\n" +
            "      number: 60\n" +
            "      frequency: 8\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 0\n" +
            "      frequency: 11\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 1\n" +
            "      frequency: 4\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 5\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 6\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 7\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 9\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 10\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 12\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -769\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -672\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -632\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -602\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -550\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -516\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -509\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -480\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -448\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -392\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -388\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -345\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -336\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 276\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 345\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 424\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 427\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 430\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 463\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 512\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 561\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 620\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 639\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 731\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 737\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 248\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 249\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 250\n" +
            "      frequency: 9\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 251\n" +
            "      frequency: 12\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 253\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 248\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 249\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 250\n" +
            "      frequency: 10\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 251\n" +
            "      frequency: 11\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 83\n" +
            "      frequency: 7\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 99\n" +
            "      frequency: 5\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 147\n" +
            "      frequency: 6\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 163\n" +
            "      frequency: 7\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -240\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -227\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -182\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -174\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -170\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -166\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -152\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -146\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -142\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -133\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -129\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -128\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -127\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -125\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -116\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -115\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -111\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -94\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -91\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -87\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -85\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -83\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -63\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -56\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -51\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -48\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -47\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -35\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -32\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -22\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -4\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -3\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 21\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 94\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 102\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -222\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -209\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -203\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -198\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -195\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -186\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -179\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -147\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -138\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -134\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -123\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -15\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -13\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -9\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -6\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 3\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 4\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 5\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 6\n" +
            "      frequency: 5\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 9\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 11\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 14\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 15\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 19\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 23\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 39\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 43\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 55\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 63\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 120\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 126\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 129\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 133\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 144\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 148\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"C\"\n" +
            "    toSequence: \"G\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "    isCalled: false\n" +
            "  }\n" +
            "  counts {\n" +
            "    matchesReference: false\n" +
            "    fromSequence: \"C\"\n" +
            "    toSequence: \"N\"\n" +
            "    genotypeCountForwardStrand: 0\n" +
            "    genotypeCountReverseStrand: 0\n" +
            "    isIndel: false\n" +
            "    qualityScoresForwardStrand {\n" +
            "      number: 0\n" +
            "      frequency: 11\n" +
            "    }\n" +
            "    qualityScoresReverseStrand {\n" +
            "      number: 0\n" +
            "      frequency: 8\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 31\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 55\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 65\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 93\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 95\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 104\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 140\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 171\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 180\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesForwardStrand {\n" +
            "      number: 206\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 23\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 26\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 37\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 55\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 134\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 161\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 187\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readIndicesReverseStrand {\n" +
            "      number: 197\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    isCalled: false\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 29\n" +
            "      frequency: 4\n" +
            "    }\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 60\n" +
            "      frequency: 6\n" +
            "    }\n" +
            "    readMappingQualityForwardStrand {\n" +
            "      number: 70\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    readMappingQualityReverseStrand {\n" +
            "      number: 29\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    readMappingQualityReverseStrand {\n" +
            "      number: 60\n" +
            "      frequency: 5\n" +
            "    }\n" +
            "    readMappingQualityReverseStrand {\n" +
            "      number: 70\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 1\n" +
            "      frequency: 10\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 2\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 3\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 4\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 6\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 7\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    numVariationsInReads {\n" +
            "      number: 9\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -597\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -545\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -472\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -461\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -452\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -437\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -373\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: -276\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 395\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 396\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 495\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 511\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 519\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 521\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 529\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 591\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 644\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    insertSizes {\n" +
            "      number: 688\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 252\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 254\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    targetAlignedLengths {\n" +
            "      number: 255\n" +
            "      frequency: 15\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 248\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 250\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    queryAlignedLengths {\n" +
            "      number: 251\n" +
            "      frequency: 15\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 83\n" +
            "      frequency: 6\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 99\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 147\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 163\n" +
            "      frequency: 7\n" +
            "    }\n" +
            "    pairFlags {\n" +
            "      number: 1123\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -208\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -188\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -180\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -164\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -155\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -146\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -138\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -133\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -123\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -117\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -105\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -94\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -84\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -81\n" +
            "      frequency: 2\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -63\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -59\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -58\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -52\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -49\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -48\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -45\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: -32\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 0\n" +
            "      frequency: 11\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 5\n" +
            "      frequency: 3\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 56\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 64\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 72\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 92\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 126\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsForwardStrand {\n" +
            "      number: 135\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -213\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: -184\n" +
            "      frequency: 1\n" +
            "    }\n" +
            "    distancesToReadVariationsReverseStrand {\n" +
            "      number: 0\n" +
            "      frequency: 8\n" +
            "    }\n" +
            "  }\n" +
            "  isTumor: true\n" +
            "  formattedCounts: \"sample: 0 counts A=0 T=0 C=25 G=0 N=0 FB=0 indels={ [indel count=19 C AGTCAGTCAGTCAGTCAGTCAGT/----AGTCAGTCAGTCAGTCAGT  45944850-45944874 filtered=false] }\\n\"\n" +
            "  isVariant: true\n" +
            "  trueGenotype: \"%s\"\n" +
            "}\n" +
            "reference_id: \"chr21\"\n" +
            "trueFrom: \"CAGTCAGTCAGTCAGTCAGTCAGT\"\n" +
            "genomicSequenceContext: \"CAGTCAGTCACTCTTTGAGGCAGTCAGTCAGTCAGTCAGTC\"\n" +
            "trueGenotype: \"%s\"";
    static String[] records = {
            createRecord("C/C----AGTCAGTCAGTCAGTCAGT"),
            createRecord("C/C----AGTCAGTCAGTCAGTCAGTAGTAGTAGTAGT")
    };

    private static String createRecord(String trueGenotype) {
        return String.format(RECORD_BASE, trueGenotype, trueGenotype);
    }
}