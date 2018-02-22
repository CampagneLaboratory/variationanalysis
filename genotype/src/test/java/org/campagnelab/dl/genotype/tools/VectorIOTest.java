package org.campagnelab.dl.genotype.tools;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.framework.tools.ExportTensorArguments;
import org.campagnelab.dl.framework.tools.ExportTensors;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV37;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.junit.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

public class VectorIOTest {

    @Test
    public void vecFileIO() throws IOException {
        try (
                SequenceBaseInformationWriter baseInformationWriter = new SequenceBaseInformationWriter("testWrite")
        ) {
            for (String record : records) {
                final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
                TextFormat.getParser().merge(record, builder);
                baseInformationWriter.appendEntry(builder.build());
            }
        } catch (IOException e) {
            throw new RuntimeException("Couldn't write sbi file");
        }
        ExportTensorArguments exportTensorArguments = new ExportTensorArguments();
        exportTensorArguments.vecFileType = "binary";
        exportTensorArguments.trainingSets = Collections.singletonList("testWrite.sbi");
        exportTensorArguments.inputNamesToExport = Collections.singleton("input");
        exportTensorArguments.outputNamesToExport = new HashSet<>(Arrays.asList("softmaxGenotype", "metaData"));
        exportTensorArguments.sampleTypes = Collections.singletonList("sampleType");
        exportTensorArguments.sampleNames = Collections.singletonList("sampleName");
        exportTensorArguments.featureMapperClassname = GenotypeMapperV37.class.getCanonicalName();
        ExportTensors exportTensors = new ExportTensorsG();
        exportTensors.arguments = exportTensorArguments;
        exportTensors.execute();


    }

    String[] records = {
            "reference_index: 19\n" +
                    "position: 43442793\n" +
                    "mutated: false\n" +
                    "referenceBase: \"C\"\n" +
                    "samples {\n" +
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
                    "    genotypeCountForwardStrand: 13\n" +
                    "    genotypeCountReverseStrand: 12\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 30\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 31\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 32\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 33\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 35\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 39\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 33\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 34\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 36\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 39\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 41\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 2\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 7\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 21\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 21\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 13\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 1\n" +
                    "      frequency: 25\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 0\n" +
                    "      frequency: 25\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 22\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 23\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 27\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 29\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 33\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 35\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 41\n" +
                    "      frequency: 36\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 22\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 23\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 27\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 33\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 35\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 41\n" +
                    "      frequency: 18\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 83\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 99\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 147\n" +
                    "      frequency: 7\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 163\n" +
                    "      frequency: 5\n" +
                    "    }\n" +
                    "    distancesToReadVariationsForwardStrand {\n" +
                    "      number: 0\n" +
                    "      frequency: 13\n" +
                    "    }\n" +
                    "    distancesToReadVariationsReverseStrand {\n" +
                    "      number: 0\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    queryPositions {\n" +
                    "      number: 0\n" +
                    "      frequency: 25\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 2\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 7\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 15\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 19\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 20\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 21\n" +
                    "      frequency: 8\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 2\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 20\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 21\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "  }\n" +
                    "  counts {\n" +
                    "    matchesReference: true\n" +
                    "    fromSequence: \"C\"\n" +
                    "    toSequence: \"C\"\n" +
                    "    genotypeCountForwardStrand: 12\n" +
                    "    genotypeCountReverseStrand: 14\n" +
                    "    isIndel: false\n" +
                    "    qualityScoresForwardStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 12\n" +
                    "    }\n" +
                    "    qualityScoresReverseStrand {\n" +
                    "      number: 40\n" +
                    "      frequency: 14\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 19\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesForwardStrand {\n" +
                    "      number: 21\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 10\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readIndicesReverseStrand {\n" +
                    "      number: 21\n" +
                    "      frequency: 13\n" +
                    "    }\n" +
                    "    isCalled: true\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 37\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityForwardStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 11\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    readMappingQualityReverseStrand {\n" +
                    "      number: 60\n" +
                    "      frequency: 13\n" +
                    "    }\n" +
                    "    numVariationsInReads {\n" +
                    "      number: 0\n" +
                    "      frequency: 26\n" +
                    "    }\n" +
                    "    insertSizes {\n" +
                    "      number: 0\n" +
                    "      frequency: 26\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 24\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 26\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 29\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 30\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 34\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 37\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 39\n" +
                    "      frequency: 4\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 40\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    targetAlignedLengths {\n" +
                    "      number: 41\n" +
                    "      frequency: 34\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 24\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 26\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 29\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 30\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 34\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 37\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 39\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 40\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    queryAlignedLengths {\n" +
                    "      number: 41\n" +
                    "      frequency: 17\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 83\n" +
                    "      frequency: 6\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 99\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 147\n" +
                    "      frequency: 7\n" +
                    "    }\n" +
                    "    pairFlags {\n" +
                    "      number: 163\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    queryPositions {\n" +
                    "      number: 0\n" +
                    "      frequency: 26\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 5\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 9\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 16\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 18\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 19\n" +
                    "      frequency: 2\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 20\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    distanceToStartOfRead {\n" +
                    "      number: 21\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 3\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 10\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 13\n" +
                    "      frequency: 1\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 20\n" +
                    "      frequency: 10\n" +
                    "    }\n" +
                    "    distanceToEndOfRead {\n" +
                    "      number: 21\n" +
                    "      frequency: 13\n" +
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
                    "    isCalled: false\n" +
                    "  }\n" +
                    "  isTumor: true\n" +
                    "  formattedCounts: \"sample: 0 counts A=0 T=25 C=26 G=0 N=0 FB=0 indels={ [] }\\n\"\n" +
                    "  isVariant: true\n" +
                    "}\n" +
                    "trueGenotype: \"C|T\"\n" +
                    "reference_id: \"chr19\"\n" +
                    "genomicSequenceContext: \"GAGTTCACGAGGGGAAAGCGCCCAGTGTGGCGTTGCGTGCA\""
    };
}


