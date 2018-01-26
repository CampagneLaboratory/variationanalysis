package org.campagnelab.dl.framework.tools;

import it.unimi.dsi.fastutil.io.FastBufferedOutputStream;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.UncheckedIOException;

public class VectorWriterBinary extends VectorWriter {
    private DataOutputStream outputStream;

    public VectorWriterBinary(String basename) throws IOException {
        super(basename);
        outputStream = new DataOutputStream(
                new FastBufferedOutputStream(new FileOutputStream(basename + ".vec"))
        );
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
        try {
            outputStream.writeInt(vectorLine.getSampleId());
            outputStream.writeLong(vectorLine.getExampleId());
            outputStream.writeInt(vectorLine.getVectorId());
            outputStream.writeInt(vectorLine.getVectorElements().size());
            for (float vectorElement : vectorLine.getVectorElements()) {
                outputStream.writeFloat(vectorElement);
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Unable to write bytes for vector line to file", e);
        }
    }
}
