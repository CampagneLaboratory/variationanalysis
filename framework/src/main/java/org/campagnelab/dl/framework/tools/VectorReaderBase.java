package org.campagnelab.dl.framework.tools;

import java.io.Closeable;
import java.io.IOException;

public abstract class VectorReaderBase implements Closeable {
    private final String inputPath;
    private final VectorWriter.VectorProperties vectorProperties;

    public VectorReaderBase(String inputPath, VectorWriter.VectorProperties vectorProperties) throws IOException {
        this.inputPath = inputPath;
        this.vectorProperties = vectorProperties;
    }

    public abstract VectorWriter.VectorLine getNextVectorLine() throws IOException;
}
