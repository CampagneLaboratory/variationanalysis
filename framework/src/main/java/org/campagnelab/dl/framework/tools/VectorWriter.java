package org.campagnelab.dl.framework.tools;

import com.google.gson.Gson;
import com.google.gson.stream.JsonWriter;
import org.nd4j.linalg.api.iter.NdIndexIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public abstract class VectorWriter implements Closeable {
    final static private int majorVersion = 0;
    final static private int minorVersion = 2;
    private JsonWriter outputFileVectorProperties;
    private int numRecords = -1;
    private List<VectorProperties.VectorPropertiesSample> sampleInfos;
    private List<VectorProperties.VectorPropertiesVector> vectorInfos;
    private Map<String, Integer> vectorNameToId;
    private Map<Integer, String> vectorIdToName;
    private Map<Integer, int[]> vectorIdToDimension;
    private int currVectorIndex = 0;
    private boolean addVectorInfoCalled = false;
    private String domainDescriptor;
    private String featureMapper;
    private String[] inputFiles;

    public VectorWriter(String basename) throws IOException {
        outputFileVectorProperties = new JsonWriter(new PrintWriter(basename + ".vecp",
                "UTF-8"));
        outputFileVectorProperties.setIndent("    ");
        sampleInfos = new LinkedList<>();
        vectorInfos = new LinkedList<>();
        vectorNameToId = new HashMap<>();
        vectorIdToName = new HashMap<>();
        vectorIdToDimension = new HashMap<>();
    }

    @Override
    public void close() throws IOException {
        Gson gson = new Gson();
        VectorProperties.VectorPropertiesSample[] sampleInfoArray = new VectorProperties.VectorPropertiesSample[
                sampleInfos.size()];
        sampleInfoArray = sampleInfos.toArray(sampleInfoArray);
        VectorProperties.VectorPropertiesVector[] vectorInfoArray;
        if (addVectorInfoCalled) {
            vectorInfoArray = new VectorProperties.VectorPropertiesVector[
                    vectorInfos.size()];
            vectorInfoArray = vectorInfos.toArray(vectorInfoArray);
        } else {
            vectorInfoArray = new VectorProperties.VectorPropertiesVector[
                    currVectorIndex];
            for (int i = 0; i < currVectorIndex; i++) {
                vectorInfoArray[i] = new VectorProperties.VectorPropertiesVector(vectorIdToName.get(i),
                        "float32", vectorIdToDimension.get(i));
            }
        }
        VectorProperties vectorProperties = new VectorProperties(getFileType(), majorVersion, minorVersion,
                numRecords, sampleInfoArray, vectorInfoArray, domainDescriptor, featureMapper, inputFiles);
        gson.toJson(vectorProperties, VectorProperties.class, outputFileVectorProperties);
        this.outputFileVectorProperties.close();
    }

    public void setNumRecords(int numRecords) {
        this.numRecords = numRecords;
    }

    public void setDomainDescriptor(String domainDescriptor) {
        this.domainDescriptor = domainDescriptor;
    }

    public void setFeatureMapper(String featureMapper) {
        this.featureMapper = featureMapper;
    }

    public void setInputFiles(String[] inputFiles) {
        this.inputFiles = inputFiles;
    }

    public void addSampleInfo(String sampleType, String sampleName) {
        this.sampleInfos.add(new VectorProperties.VectorPropertiesSample(sampleType, sampleName));
    }

    public void addVectorInfo(String vectorName, String vectorType, int[] vectorDimension) {
        addVectorInfoCalled = true;
        vectorNameToId.put(vectorName, vectorInfos.size());
        this.vectorInfos.add(new VectorProperties.VectorPropertiesVector(vectorName, vectorType, vectorDimension));
    }

    public void appendMdsList(List<MultiDataSet> multiDataSetList, int[] inputIndices, int[] outputIndices,
                              String[] inputNames, String[] outputNames, long startExampleIndex) {
        int numExamplesInBatch = multiDataSetList.get(0).getFeatures(inputIndices[0]).rows();
        if (startExampleIndex + numExamplesInBatch > numRecords) {
            throw new IllegalArgumentException("Example ID exceeds number of records");
        }
        for (int currExampleInBatch = 0; currExampleInBatch < numExamplesInBatch; currExampleInBatch++) {
            int sampleMdsIndex = 0;
            for (MultiDataSet multiDataSetAtSample : multiDataSetList) {
                writeLinesForExample(multiDataSetAtSample, inputIndices, inputNames, sampleMdsIndex,
                        startExampleIndex, currExampleInBatch, numExamplesInBatch, true);
                writeLinesForExample(multiDataSetAtSample, outputIndices, outputNames, sampleMdsIndex,
                        startExampleIndex, currExampleInBatch, numExamplesInBatch, false);
                sampleMdsIndex++;
            }
        }
    }

    private int[] getShape(INDArray values) {
        // Store row and column vectors as 1-dimensional
        if (values.isRowVector()) {
            return new int[]{values.columns()};
        } else if (values.isColumnVector()) {
            return new int[]{values.rows()};
        } else {
            return values.shape();
        }
    }

    private void writeLinesForExample(MultiDataSet multiDataSet, int[] indices, String[] names,
                                      int sampleIndex, long startExampleIndex,
                                      int currExampleIndexInBatch, int numExamplesInBatch,
                                      boolean isForFeatures) {
        for (int index : indices) {
            INDArray allValuesAtIndex = isForFeatures
                    ? multiDataSet.getFeatures(index)
                    : multiDataSet.getLabels(index);
            if (allValuesAtIndex.rows() != numExamplesInBatch) {
                throw new RuntimeException("Mismatched mds dimensions for batch size");
            }
            INDArray currExampleValuesAtIndex = allValuesAtIndex.getRow(currExampleIndexInBatch);
            String vectorName = names[index];
            int vectorId;
            if (vectorNameToId.get(vectorName) == null) {
                vectorNameToId.put(vectorName, currVectorIndex);
                vectorIdToName.put(currVectorIndex, vectorName);
                vectorIdToDimension.put(currVectorIndex, getShape(currExampleValuesAtIndex));

                vectorId = currVectorIndex++;
            } else {
                vectorId = vectorNameToId.get(vectorName);
                String vectorCachedName = vectorIdToName.get(vectorId);
                if (!vectorName.equals(vectorCachedName)) {
                    throw new RuntimeException(String.format("Vector name mismatch for vector id %d", vectorId));
                }
                int[] vectorCachedDimensions = vectorIdToDimension.get(vectorId);
                if (!Arrays.equals(getShape(currExampleValuesAtIndex), vectorCachedDimensions)) {
                    throw new RuntimeException(String.format("Vector dimension mismatch for vector id %d", vectorId));
                }
            }
            writeVectorLine(new VectorLine(sampleIndex, startExampleIndex + currExampleIndexInBatch, vectorId,
                        getVectorElementsFromArray(currExampleValuesAtIndex)));
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
        private long exampleId;
        private int vectorId;
        private List<Float> vectorElements;

        public VectorLine(int sampleId, long exampleId, int vectorId, List<Float> vectorElements) {
            this.sampleId = sampleId;
            this.exampleId = exampleId;
            this.vectorId = vectorId;
            this.vectorElements = vectorElements;
        }

        public int getSampleId() {
            return sampleId;
        }

        public long getExampleId() {
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
        private int majorVersion;
        private int minorVersion;
        private String fileType;
        private long numRecords;
        private String domainDescriptor;
        private String featureMapper;
        private String[] inputFiles;
        private VectorPropertiesSample[] samples;
        private VectorPropertiesVector[] vectors;

        VectorProperties(String fileType, int majorVersion, int minorVersion, long numRecords,
                         VectorPropertiesSample[] samples, VectorPropertiesVector[] vectors,
                         String domainDescriptor, String featureMapper, String[] inputFiles) {
            this.fileType = fileType;
            this.majorVersion = majorVersion;
            this.minorVersion = minorVersion;
            this.samples = samples;
            this.vectors = vectors;
            this.numRecords = numRecords;
            this.domainDescriptor = domainDescriptor;
            this.featureMapper = featureMapper;
            this.inputFiles = inputFiles;
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

        public long getNumRecords() { return numRecords; }

        public String getDomainDescriptor() {
            return domainDescriptor;
        }

        public String getFeatureMapper() {
            return featureMapper;
        }

        public String[] getInputFiles() {
            return inputFiles;
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
