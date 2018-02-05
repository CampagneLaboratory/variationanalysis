package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

public class PyTorchModelClient {
    private final String modelCommand;
    private final String modelName;
    private final String checkpointDir;
    private final int sampleId;
    private final String[] vectorOutputNames;
    private final String sampleType;
    private final String sampleName;
    private final String[] vectorInputNames;
    private int batchId;

    enum ModelType {
        GENOTYPE,
        SOMATIC
    }

    public PyTorchModelClient(String modelName, String checkpointDir, ModelType modelType,
                              DomainDescriptor domainDescriptor) {
        this(modelName, checkpointDir, modelType, domainDescriptor, 0);
    }

    public PyTorchModelClient(String modelName, String checkpointDir, ModelType modelType,
                              DomainDescriptor domainDescriptor, int sampleId) {
        this.modelName = modelName;
        this.checkpointDir = checkpointDir;
        switch (modelType) {
            case SOMATIC:
                modelCommand = "infer_somatic";
                sampleType = "somaticType";
                sampleName = "somaticName";
                break;
            case GENOTYPE:
                modelCommand = "infer_genotype";
                sampleType = "genotypeType";
                sampleName = "genotypeName";
                break;
            default:
                throw new UnsupportedOperationException();
        }
        vectorOutputNames = domainDescriptor.getComputationalGraph().getOutputNames();
        vectorInputNames = domainDescriptor.getComputationalGraph().getInputNames();
        this.sampleId = sampleId;
        batchId = 0;
    }

    public INDArray[] predict(MultiDataSet multiDataSet, int miniBatchSize) throws IOException {
        // Create temporary file for mini batch to use for inference in genotypetensors
        File tempFile = File.createTempFile("mds" + batchId++, "vec");
        tempFile.deleteOnExit();
        String tempPath = tempFile.getCanonicalPath();
        String tempBaseName = FilenameUtils.getBaseName(tempPath);
        // Cache minibatch features to vector file
        try (
                VectorWriter vectorWriter = new VectorWriterBinary(tempBaseName)
        ) {
            vectorWriter.addSampleInfo(sampleType, sampleName);
            int batchSize = multiDataSet.getFeatures()[0].rows();
            for (int i = 0; i < batchSize; i++) {
                for (int j = 0; j < vectorInputNames.length; j++) {
                    FloatArrayList inputs = new FloatArrayList(multiDataSet.getFeatures(j).ravel('c').data().asFloat());
                    vectorWriter.writeVectorLine(new VectorWriter.VectorLine(sampleId, i, j, inputs));
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to map vectors for PyTorch inference");
        }
        // Call inference for cached file
        ProcessBuilder processBuilder = new ProcessBuilder("zerorpc", "tcp://0:1234",
                modelCommand, modelName, checkpointDir, tempPath);
        processBuilder.redirectErrorStream(true);
        Process process = processBuilder.start();
        BufferedReader processOutputReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line;
        int linesOutputted = 0;
        // Get back name of output vector file for prediction
        String outputFileName = "";
        while ((line = processOutputReader.readLine()) != null) {
            if (linesOutputted == 0) {
                outputFileName = line;
            }
            linesOutputted++;
        }
        VectorReader vectorReader = new VectorReader(outputFileName, sampleId, vectorOutputNames);
        VectorReader.RecordVectors recordVectors = vectorReader.getNextBatch(miniBatchSize);
        return recordVectors.getVectors();
    }


}
