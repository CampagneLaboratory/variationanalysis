package org.campagnelab.dl.genotype.tools;

import com.google.gson.Gson;
import com.google.gson.stream.JsonWriter;
import org.nd4j.linalg.api.iter.NdIndexIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.MultiDataSet;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public abstract class VectorWriter implements Closeable {
    private JsonWriter outputFileVectorProperties;
    private int majorVersion = -1;
    private int minorVersion = -1;
    private List<VectorProperties.VectorPropertiesSample> sampleInfos;
    private List<VectorProperties.VectorPropertiesVector> vectorInfos;
    private Map<String, Integer> vectorNameToId;

    public VectorWriter(String basename) throws IOException {
        outputFileVectorProperties = new JsonWriter(new PrintWriter(basename + ".vecp",
                "UTF-8"));
        outputFileVectorProperties.setIndent("    ");
        sampleInfos = new LinkedList<>();
        vectorInfos = new LinkedList<>();
        vectorNameToId = new HashMap<>();
    }

    @Override
    public void close() throws IOException {
        Gson gson = new Gson();
        if ((majorVersion == -1) || (minorVersion == -1)) {
            majorVersion = 0;
            minorVersion = 1;
        }
        VectorProperties.VectorPropertiesSample[] sampleInfoArray = new VectorProperties.VectorPropertiesSample[
                sampleInfos.size()];
        VectorProperties.VectorPropertiesVector[] vectorInfoArray = new VectorProperties.VectorPropertiesVector[
                vectorInfos.size()];
        sampleInfoArray = sampleInfos.toArray(sampleInfoArray);
        vectorInfoArray = vectorInfos.toArray(vectorInfoArray);
        VectorProperties vectorProperties = new VectorProperties(getFileType(), majorVersion, minorVersion,
                sampleInfoArray, vectorInfoArray);
        gson.toJson(vectorProperties, VectorProperties.class, outputFileVectorProperties);
        this.outputFileVectorProperties.close();
    }

    public void setSpecVersionNumber(int majorVersion, int minorVersion) {
        this.majorVersion = majorVersion;
        this.minorVersion = minorVersion;
    }

    public void addSampleInfo(String sampleType, String sampleName) {
        this.sampleInfos.add(new VectorProperties.VectorPropertiesSample(sampleType, sampleName));
    }

    public void addVectorInfo(String vectorName, String vectorType, int[] vectorDimension) {
        vectorNameToId.put(vectorName, vectorInfos.size());
        this.vectorInfos.add(new VectorProperties.VectorPropertiesVector(vectorName, vectorType, vectorDimension));
    }

    public void appendMds(MultiDataSet multiDataSet, int[] inputIndices, int[] outputIndices, String[] inputNames,
                          String[] outputNames) {
        if ((inputIndices.length != inputNames.length) || (outputIndices.length != outputNames.length)) {
            throw new IllegalArgumentException("Lengths of input indices and names, and output indices and names," +
                    " must be equal");
        }
        for (int i = 0; i < inputIndices.length; i++) {
            writeVectorLine(new VectorLine(0, 0, vectorNameToId.get(inputNames[i]),
                    getVectorElementsFromArray(multiDataSet.getFeatures(inputIndices[i]))));
        }
        for (int i = 0; i < outputIndices.length; i++) {
            writeVectorLine(new VectorLine(0, 0, vectorNameToId.get(outputNames[i]),
                    getVectorElementsFromArray(multiDataSet.getFeatures(outputIndices[i]))));
        }

    }

    private List<Float> getVectorElementsFromArray(INDArray vectorArray) {
        List<Float> vectorElements = new LinkedList<>();
        NdIndexIterator ndIndexIterator = new NdIndexIterator('c', vectorArray.shape());
        while (ndIndexIterator.hasNext()) {
            vectorElements.add(vectorArray.getFloat(ndIndexIterator.next()));
        }
        return vectorElements;
    }

    public abstract String getFileType();

    public abstract void writeVectorLine(VectorLine vectorLine);

    static class VectorLine {
        private int sampleId;
        private int exampleId;
        private int vectorId;
        private List<Float> vectorElements;

        public VectorLine(int sampleId, int exampleId, int vectorId, List<Float> vectorElements) {
            this.sampleId = sampleId;
            this.exampleId = exampleId;
            this.vectorId = vectorId;
            this.vectorElements = vectorElements;
        }

        public int getSampleId() {
            return sampleId;
        }

        public int getExampleId() {
            return exampleId;
        }

        public int getVectorId() {
            return vectorId;
        }

        public List<Float> getVectorElements() {
            return vectorElements;
        }
    }

    static class VectorProperties {
        private String fileType;
        private int majorVersion;
        private int minorVersion;
        private VectorPropertiesSample[] samples;
        private VectorPropertiesVector[] vectors;

        VectorProperties(String fileType, int majorVersion, int minorVersion, VectorPropertiesSample[] samples,
                         VectorPropertiesVector[] vectors) {
            this.fileType = fileType;
            this.majorVersion = majorVersion;
            this.minorVersion = minorVersion;
            this.samples = samples;
            this.vectors = vectors;
        }

        public String getFileType() {
            return fileType;
        }

        public int getMajorVersion() {
            return majorVersion;
        }

        public int getMinorVersion() {
            return minorVersion;
        }

        public VectorPropertiesSample[] getSamples() {
            return samples;
        }

        public VectorPropertiesVector[] getVectors() {
            return vectors;
        }

        static class VectorPropertiesSample {
            private String sampleType;
            private String sampleName;

            VectorPropertiesSample(String sampleType, String sampleName) {
                this.sampleType = sampleType;
                this.sampleName = sampleName;
            }

            public String getSampleType() {
                return sampleType;
            }

            public String getSampleName() {
                return sampleName;
            }
        }

        static class VectorPropertiesVector {
            private String vectorName;
            private String vectorType;
            private int[] vectorDimension;

            VectorPropertiesVector(String vectorName, String vectorType, int[] vectorDimension) {
                this.vectorName = vectorName;
                this.vectorType = vectorType;
                this.vectorDimension = vectorDimension;
            }

            public String getVectorName() {
                return vectorName;
            }

            public String getVectorType() {
                return vectorType;
            }

            public int[] getVectorDimension() {
                return vectorDimension;
            }
        }
    }
}
