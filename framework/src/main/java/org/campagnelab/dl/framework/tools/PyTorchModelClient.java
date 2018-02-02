package org.campagnelab.dl.framework.tools;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
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
    private int exampleId;

    enum ModelType {
        GENOTYPE,
        SOMATIC
    }

    public PyTorchModelClient(String modelName, String checkpointDir, ModelType modelType) {
        this(modelName, checkpointDir, modelType, 0);
    }

    public PyTorchModelClient(String modelName, String checkpointDir, ModelType modelType, int sampleId) {
        this.modelName = modelName;
        this.checkpointDir = checkpointDir;
        switch (modelType) {
            case SOMATIC:
                modelCommand = "infer_somatic";
                vectorOutputNames = new String[]{"isBaseMutated", "somaticFrequency"};
                sampleType = "somaticType";
                sampleName = "somaticName";
                break;
            case GENOTYPE:
                modelCommand = "infer_genotype";
                vectorOutputNames = new String[]{"softmaxGenotype"};
                sampleType = "genotypeType";
                sampleName = "genotypeName";
                break;
            default:
                throw new UnsupportedOperationException();
        }
        this.sampleId = sampleId;
        exampleId = 0;
    }

    public INDArray[] predict(MultiDataSet multiDataSet) throws IOException {
        File tempFile = File.createTempFile("mds" + exampleId++, "vec");
        tempFile.deleteOnExit();
        String tempPath = tempFile.getCanonicalPath();
        String tempBaseName = FilenameUtils.getBaseName(tempPath);
        try (
                VectorWriter vectorWriter = new VectorWriterBinary(tempBaseName)
        ) {
            vectorWriter.addSampleInfo(sampleType, sampleName);
        } catch (IOException e) {
            throw new RuntimeException("Unable to map vectors for PyTorch inference");
        }

        ProcessBuilder processBuilder = new ProcessBuilder("zerorpc", "tcp://0:1234",
                modelCommand, modelName, checkpointDir, tempPath);
        processBuilder.redirectErrorStream(true);
        Process process = processBuilder.start();
        BufferedReader processOutputReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        String line;
        int linesOutputted = 0;
        String outputFileName = "";
        while ((line = processOutputReader.readLine()) != null) {
            if (linesOutputted == 0) {
                outputFileName = line;
            }
            linesOutputted++;
        }
        VectorReader vectorReader = new VectorReader(outputFileName, sampleId, vectorOutputNames);
        VectorReader.RecordVectors recordVectors;
        return null;
    }


}
