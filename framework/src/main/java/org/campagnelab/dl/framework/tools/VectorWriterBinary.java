package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public class VectorWriterBinary extends VectorWriter {
    // Size in bytes of different fields needed for binary format
    private static final int EXAMPLE_ID_SIZE = 8;
    private static final int VECTOR_ID_SIZE = 4;
    private static final int SAMPLE_ID_SIZE = 4;
    private static final int VECTOR_LENGTH_SIZE = 4;
    private static final int HEADER_SIZE = EXAMPLE_ID_SIZE + VECTOR_ID_SIZE + SAMPLE_ID_SIZE + VECTOR_LENGTH_SIZE;
    private static final int VECTOR_ELEMENT_SIZE = 4;

    private FastBufferedOutputStream outputStream;

    public VectorWriterBinary(String basename) throws IOException {
        super(basename);
        outputStream = new FastBufferedOutputStream(new FileOutputStream(basename + ".vec"));
    }

    @Override
    public void close() throws IOException {
        super.close();
        outputStream.close();
    }

    @Override
    public String getFileType() {
        return "binary";
    }

    @Override
    public void writeVectorLine(VectorLine vectorLine) {
        final int bytesNeeded = HEADER_SIZE + (vectorLine.getVectorElements().size() * VECTOR_ELEMENT_SIZE);
        ByteBuffer vectorLineByteBuffer = ByteBuffer.allocate(bytesNeeded);
        vectorLineByteBuffer.order(ByteOrder.BIG_ENDIAN);
        vectorLineByteBuffer.putInt(vectorLine.getSampleId());
        vectorLineByteBuffer.putLong(vectorLine.getExampleId());
        vectorLineByteBuffer.putInt(vectorLine.getVectorId());
        vectorLineByteBuffer.putInt(vectorLine.getVectorElements().size());
        for (float vectorElement : vectorLine.getVectorElements()) {
            vectorLineByteBuffer.putFloat(vectorElement);
        }
        try {
            outputStream.write(vectorLineByteBuffer.array());
        } catch (IOException e) {
            throw new UncheckedIOException("Unable to write bytes for vector line to file", e);
        }
    }
}
