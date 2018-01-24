package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.tools.VectorWriter;
import org.campagnelab.goby.compression.GzipOutputStreamWithCustomLevel;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

public class VectorWriterText extends VectorWriter {
    private PrintWriter outputFileVector;

    public VectorWriterText(String basename) throws IOException {
        super(basename);
        outputFileVector = new PrintWriter(new GzipOutputStreamWithCustomLevel(1, new FileOutputStream(basename+".vec")));
    }

    @Override
    public void close() throws IOException {
        super.close();
        outputFileVector.close();
    }

    @Override
    public String getFileType() {
        return "gzipped+text";
    }

    @Override
    public void writeVectorLine(VectorLine vectorLine) {
        outputFileVector.append(Integer.toString(vectorLine.getSampleId()))
                .append(" ")
                .append(Long.toString(vectorLine.getExampleId()))
                .append(" ")
                .append(Integer.toString(vectorLine.getVectorId()));
        for (float vectorValue : vectorLine.getVectorElements()) {
            outputFileVector.append(" ")
                    .append(Float.toString(vectorValue));
        }
        outputFileVector.append("\n");
    }
}
