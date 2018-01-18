package org.campagnelab.dl.framework.tools;

import org.campagnelab.dl.framework.tools.VectorWriter;

import java.io.IOException;
import java.io.PrintWriter;

public class VectorWriterText extends VectorWriter {
    private PrintWriter outputFileVector;

    public VectorWriterText(String basename) throws IOException {
        super(basename);
        outputFileVector = new PrintWriter(basename + ".vec", "UTF-8");
    }

    @Override
    public void close() throws IOException {
        super.close();
        outputFileVector.close();
    }

    @Override
    public String getFileType() {
        return "text";
    }

    @Override
    public void writeVectorLine(VectorLine vectorLine) {
        outputFileVector.append(Integer.toString(vectorLine.getSampleId()))
                .append(" ")
                .append(Integer.toString(vectorLine.getExampleId()))
                .append(" ")
                .append(Integer.toString(vectorLine.getVectorId()));
        for (float vectorValue : vectorLine.getVectorElements()) {
            outputFileVector.append(" ")
                    .append(Float.toString(vectorValue));
        }
        outputFileVector.append("\n");
    }
}
