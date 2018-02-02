package org.campagnelab.dl.framework.tools;

import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;
import it.unimi.dsi.fastutil.longs.LongArrayList;
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
import java.util.*;

public class VectorReader implements Closeable {
    private final VectorWriter.VectorProperties vectorProperties;
    private final VectorReaderBase delegateReader;
    private final Set<Pair<Integer, Integer>> sampleVectorIds;
    private final Set<Long> processedExampleIds;
    private final JsonReader propertiesReader;
    private final int sampleId;
    private final Map<Integer, Integer> vectorIds;
    private final Map<Integer, int[]> vectorDimensions;
    private final String[] vectorNames;
    private final boolean returnExampleIds;
    private final int[] vectorIdArray;

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
        vectorIds = new HashMap<>();
        vectorDimensions = new HashMap<>();
        vectorIdArray = new int[vectorNames.length];
        int vectorIdArrayIdx = 0;
        for (String vectorName : vectorNames) {
            int i = 0;
            for (VectorWriter.VectorProperties.VectorPropertiesVector vectorInfo : vectorProperties.getVectors()) {
                if (vectorInfo.getVectorName().equals(vectorName)) {
                    vectorIds.put(i, vectorIdArrayIdx);
                    vectorIdArray[vectorIdArrayIdx++] = i;
                    vectorDimensions.put(i, vectorInfo.getVectorDimension());
                    break;
                }
                i++;
            }
        }
        if (vectorIds.size() != vectorNames.length) {
            throw new RuntimeException("Vector names not found in vector properties file");
        }
        this.vectorNames = vectorNames;
        this.returnExampleIds = returnExampleIds;
    }

    public RecordVectors getNextBatch() {
        return this.getNextBatch(1);
    }

    public RecordVectors getNextBatch(int batchSize) {
        LongArrayList nextExampleIds = new LongArrayList(batchSize);
        ObjectArrayList<INDArray[]> nextExampleArrays = new ObjectArrayList<>(batchSize);
        Pair<Long, INDArray[]> nextExamplePair;
        for (int i = 0; i < batchSize; i++) {
            if ((nextExamplePair = getNextExampleList()) != null) {
                nextExampleIds.add(nextExamplePair.getLeft());
                nextExampleArrays.add(nextExamplePair.getRight());
            } else {
                break;
            }
        }
        nextExampleIds.trim();
        nextExampleArrays.trim();
        INDArray[] indArrays = new INDArray[vectorIds.size()];
        for (int i = 0; i < indArrays.length; i++) {
            int[] dims = vectorDimensions.get(vectorIdArray[i]);
            int[] shape = new int[dims.length + 1];
            shape[0] = nextExampleIds.size();
            System.arraycopy(shape, 0, dims, 1, shape.length);
            indArrays[i] = Nd4j.create(shape, 'c');
        }
        long[] nextExampleIdArray = new long[nextExampleIds.size()];
        nextExampleIdArray = nextExampleIds.toArray(nextExampleIdArray);
        for (int i = 0; i < nextExampleArrays.size(); i++) {
            INDArray[] nextExampleINDArray = nextExampleArrays.get(i);
            for (int j = 0; j < nextExampleINDArray.length; j++) {
                indArrays[j].putRow(i, nextExampleINDArray[j]);
            }
        }
        return new RecordVectors(nextExampleIdArray, vectorNames, indArrays);
    }

    private Pair<Long, INDArray[]> getNextExampleList() {
        INDArray[] vectorArrays = new INDArray[vectorIds.size()];
        Set<Pair<Integer, Integer>> processedVectorSampleIds = new HashSet<>();
        long currExampleId = -1;
        int filledSlots = 0;
        try {
            for (int i = 0; i < sampleVectorIds.size(); i++) {
                VectorWriter.VectorLine vectorLine = delegateReader.getNextVectorLine();
                if (currExampleId == -1) currExampleId = vectorLine.getExampleId();
                if (processedExampleIds != null) {
                    if (processedExampleIds.contains(currExampleId)) {
                        throw new RuntimeException(String.format("Example ID %d already processed", currExampleId));
                    }
                }
                Integer vectorIndexInArray = vectorIds.get(vectorLine.getVectorId());
                if ((vectorIndexInArray != null) && vectorLine.getSampleId() == sampleId) {
                    float[] vectorElementsArray = new float[vectorLine.getVectorElements().size()];
                    vectorElementsArray = vectorLine.getVectorElements().toArray(vectorElementsArray);
                    int[] vectorShape = vectorProperties.getVectors()[vectorLine.getVectorId()].getVectorDimension();
                    vectorArrays[vectorIndexInArray] = Nd4j.create(vectorElementsArray, vectorShape, 'c');
                    filledSlots++;
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
            if (!(filledSlots == vectorIds.size())) {
                throw new RuntimeException("Missing vectors");
            }
            return new ImmutablePair<>(currExampleId, vectorArrays);
        } catch (IOException e) {
            return null;
        }
    }

    @Override
    public void close() throws IOException {
        delegateReader.close();
        propertiesReader.close();
    }

    public static class RecordVectors {
        private long[] exampleIds;
        private String[] vectorNames;
        private INDArray[] vectors;

        public RecordVectors(String[] vectorNames, INDArray[] vectors) {
            this(new long[]{0L}, vectorNames, vectors);
        }

        public RecordVectors(long[] exampleIds, String[] vectorNames, INDArray[] vectors) {
            this.exampleIds = exampleIds;
            this.vectorNames = vectorNames;
            this.vectors = vectors;
        }

        public boolean hasExampleIds() {
            return exampleIds != null;
        }

        public long[] getExampleIds() {
            if (!hasExampleIds()) throw new UnsupportedOperationException("No example ID present");
            return exampleIds;
        }

        public String[] getVectorNames() {
            return vectorNames;
        }

        public INDArray[] getVectors() {
            return vectors;
        }

        @Override
        public String toString() {
            return String.format("Example IDs: %d\nVector names: %s\nVectors: %s\n",
                    Arrays.toString(exampleIds), Arrays.toString(vectorNames), Arrays.toString(vectors));
        }
    }
}
