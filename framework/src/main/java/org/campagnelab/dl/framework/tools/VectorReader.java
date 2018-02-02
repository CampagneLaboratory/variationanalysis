package org.campagnelab.dl.framework.tools;

import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.Closeable;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

public class VectorReader implements Closeable {
    private final VectorWriter.VectorProperties vectorProperties;
    private final VectorReaderBase delegateReader;
    private final Set<Pair<Integer, Integer>> sampleVectorIds;
    private final Set<Long> processedExampleIds;
    private final JsonReader propertiesReader;
    private final int sampleId;
    private final Set<Integer> vectorIds;
    private final String[] vectorNames;
    private final boolean returnExampleIds;

    public VectorWriter.VectorProperties getVectorProperties() {
        return vectorProperties;
    }

    public VectorReader(String inputPath, int sampleId, String[] vectorNames) throws IOException {
        this(inputPath, sampleId, vectorNames,true);
    }

    public VectorReader(String inputPath, int sampleId, String[] vectorNames,
                        boolean returnExampleIds) throws IOException {
        this(inputPath, sampleId, vectorNames, false, returnExampleIds);
    }

    public VectorReader(String inputPath, int sampleId, String[] vectorNames,
                        boolean assertExampleIds, boolean returnExampleIds) throws IOException {
        String propertiesFileName = FilenameUtils.removeExtension(inputPath) + ".vecp";
        propertiesReader = new JsonReader(new InputStreamReader(new FileInputStream(propertiesFileName),
                "UTF-8"));
        vectorProperties = new Gson().fromJson(propertiesReader, VectorWriter.VectorProperties.class);
        switch (vectorProperties.getFileType()) {
            case "text":
            case "gzipped+text":
                throw new UnsupportedOperationException();
            case "binary":
                delegateReader = new VectorReaderBinary(inputPath, vectorProperties);
                break;
            default:
                throw new UnsupportedOperationException();
        }
        sampleVectorIds = new HashSet<>();
        for (int i = 0; i < vectorProperties.getSamples().length; i++) {
            for (int j = 0; j < vectorProperties.getVectors().length; j++) {
                sampleVectorIds.add(new ImmutablePair<>(i, j));
            }
        }
        processedExampleIds = assertExampleIds ? new HashSet<>() : null;
        this.sampleId = sampleId;
        vectorIds = new HashSet<>();
        for (String vectorName : vectorNames) {
            int i = 0;
            for (VectorWriter.VectorProperties.VectorPropertiesVector vectorInfo : vectorProperties.getVectors()) {
                if (vectorInfo.getVectorName().equals(vectorName)) {
                    vectorIds.add(i);
                }
                i++;
            }
        }
        this.vectorNames = vectorNames;
        this.returnExampleIds = returnExampleIds;
    }

    public RecordVectors getNextExample() {
        ObjectArrayList<INDArray> vectorArrays = new ObjectArrayList<>(vectorIds.size());
        Set<Pair<Integer, Integer>> processedVectorSampleIds = new HashSet<>();
        long currExampleId = -1;
        try {
            for (int i = 0; i < sampleVectorIds.size(); i++) {
                VectorWriter.VectorLine vectorLine = delegateReader.getNextVectorLine();
                if (currExampleId == -1) currExampleId = vectorLine.getExampleId();
                if (processedExampleIds != null) {
                    if (processedExampleIds.contains(currExampleId)) {
                        throw new RuntimeException(String.format("Example ID %d already processed", currExampleId));
                    }
                }
                if (vectorIds.contains(vectorLine.getVectorId()) && vectorLine.getSampleId() == sampleId) {
                    float[] vectorElementsArray = new float[vectorLine.getVectorElements().size()];
                    vectorElementsArray = vectorLine.getVectorElements().toArray(vectorElementsArray);
                    int[] vectorShape = vectorProperties.getVectors()[vectorLine.getVectorId()].getVectorDimension();
                    vectorArrays.add(Nd4j.create(vectorElementsArray, vectorShape, 'c'));
                }
                processedVectorSampleIds.add(new ImmutablePair<>(vectorLine.getSampleId(), vectorLine.getVectorId()));
            }
            if (!sampleVectorIds.equals(processedVectorSampleIds)) {
                sampleVectorIds.removeAll(processedVectorSampleIds);
                throw new RuntimeException(
                        String.format("Missing vector index-sample index pairs for example %d: %s",
                                currExampleId, sampleVectorIds.toString())
                );
            }
            INDArray[] returnVectorArrays = new INDArray[vectorArrays.size()];
            returnVectorArrays = vectorArrays.toArray(returnVectorArrays);
            return returnExampleIds
                    ? new RecordVectors(currExampleId, vectorNames, returnVectorArrays)
                    : new RecordVectors(vectorNames, returnVectorArrays);
        } catch (IOException e) {
            return null;
        }
    }

    @Override
    public void close() throws IOException {
        delegateReader.close();
        propertiesReader.close();
    }

    public RecordVectors getNextBatch(int batchSize) {
        if (batchSize==1) {
            return getNextExample();
        }
        //TODO implement this method when batchSize larger than 1.
        return null;
    }

    public static class RecordVectors {
        private long exampleId;
        private String[] vectorNames;
        private INDArray[] vectors;

        public RecordVectors(String[] vectorNames, INDArray[] vectors) {
            this(-1, vectorNames, vectors);
        }

        public RecordVectors(long exampleId, String[] vectorNames, INDArray[] vectors) {
            this.exampleId = exampleId;
            this.vectorNames = vectorNames;
            this.vectors = vectors;
        }

        public boolean hasExampleId() {
            return exampleId != -1;
        }

        public long getExampleId() {
            if (!hasExampleId()) throw new UnsupportedOperationException("No example ID present");
            return exampleId;
        }

        public String[] getVectorNames() {
            return vectorNames;
        }

        public INDArray[] getVectors() {
            return vectors;
        }

        @Override
        public String toString() {
            return String.format("Example ID: %d\nVector names: %s\nVectors: %s\n",
                    exampleId, Arrays.toString(vectorNames), Arrays.toString(vectors));
        }
    }
}
