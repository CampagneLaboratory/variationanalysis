package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.tools.arguments.AbstractTool;

import java.io.IOException;

public class ImportTensors extends AbstractTool<ImportTensorArguments> {
    private VectorWriter.VectorProperties vectorProperties;

    public static void main(String[] args) {
        ImportTensors importTensors = new ImportTensors();
        importTensors.parseArguments(args, "ImportTensors", importTensors.createArguments());
        importTensors.execute();
    }

    @Override
    public ImportTensorArguments args() {
        return (ImportTensorArguments) arguments;
    }

    @Override
    public ImportTensorArguments createArguments() {
        return new ImportTensorArguments();
    }

    @Override
    public void execute() {
        String[] vectorNames = new String[args().vectorNames.size()];
        vectorNames = args().vectorNames.toArray(vectorNames);
        try (
                VectorReader vectorReader = new VectorReader(args().inputPath, args().sampleId, vectorNames)
        ) {
            vectorProperties = vectorReader.getVectorProperties();
            VectorReader.RecordVectors recordVectors;
            long examplesProcessed = 0;
            while ((examplesProcessed < args().importN)
                    && ((recordVectors = vectorReader.getNextBatch(args().miniBatchSize)) != null)) {
                args().processVectors().accept(recordVectors);
                examplesProcessed++;
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to read in vectors. ", e);
        }
    }

    public VectorWriter.VectorProperties getVectorProperties() {
        return vectorProperties;
    }
}
