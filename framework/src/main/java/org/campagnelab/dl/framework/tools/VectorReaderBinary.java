package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;

public class VectorReaderBinary extends VectorReaderBase {
    private DataInputStream inputStream;

    public VectorReaderBinary(String inputPath, VectorWriter.VectorProperties vectorProperties) throws IOException {
        super(inputPath, vectorProperties);
        inputStream = new DataInputStream(new FastBufferedInputStream(new FileInputStream(inputPath)));
    }

    @Override
    public VectorWriter.VectorLine getNextVectorLine() throws IOException {
        int sampleId = inputStream.readInt();
        long exampleId = inputStream.readLong();
        int vectorId = inputStream.readInt();
        int numElements = inputStream.readInt();
        FloatArrayList vectorElements = new FloatArrayList(numElements);
        for (int i = 0; i < numElements; i++) {
            vectorElements.add(inputStream.readFloat());
        }
        return new VectorWriter.VectorLine(sampleId, exampleId, vectorId, vectorElements);
    }

    @Override
    public void close() throws IOException {
        inputStream.close();
    }
}
