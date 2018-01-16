package org.campagnelab.dl.genotype.tools;

import com.google.gson.Gson;
import com.google.gson.stream.JsonWriter;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.mappers.ConfigurableFeatureMapper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.learning.domains.GenotypeDomainDescriptor;
import org.campagnelab.dl.genotype.mappers.CombinedLabelsMapper;
import org.campagnelab.dl.genotype.mappers.GenotypeMapperV37;
import org.campagnelab.dl.somatic.learning.TrainSomaticModel;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.nd4j.linalg.api.iter.NdIndexIterator;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.DataSet;
import org.campagnelab.dl.somatic.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.Properties;


public class WriteVectors extends AbstractTool<WriteVectorsArguments> {
    static private Logger LOG = LoggerFactory.getLogger(WriteVectors.class);

    public static void main(String[] args) {
        WriteVectors tool = new WriteVectors();
        tool.parseArguments(args, "WriteVectors", tool.createArguments());
        tool.execute();
    }

    @Override
    public WriteVectorsArguments createArguments() {
        return new WriteVectorsArguments();
    }

    @Override
    public void execute() {
        FeatureMapper<BaseInformationRecords.BaseInformationOrBuilder> featureMapper = new GenotypeMapperV37();
        ConfigurableFeatureMapper confMapper = (ConfigurableFeatureMapper) featureMapper;
        Properties sbiProperties;
        try {
            sbiProperties = TrainSomaticModel.getReaderProperties(args().trainingSet);
        } catch (IOException e) {
            throw new UncheckedIOException("Couldn't read in training set properties", e);
        }
        confMapper.configure(sbiProperties);
        LabelMapper<BaseInformationRecords.BaseInformation> labelMapper = new CombinedLabelsMapper();
        BaseInformationIterator baseInformationIterator;
        try {
            baseInformationIterator = new BaseInformationIterator(args().trainingSet, 1, featureMapper,
                    labelMapper);
        } catch (IOException e) {
            throw new UncheckedIOException("Unable to load input files", e);
        }
        final String outputFilenameVector = args().outputBasename + ".vec";
        final String outputFilenameVectorProperties = args().outputBasename + ".vecp";
        PrintWriter outputFileVector;
        JsonWriter outputFileVectorProperties;
        try {
            outputFileVector = new PrintWriter(outputFilenameVector, "UTF-8");
            outputFileVectorProperties = new JsonWriter(new PrintWriter(outputFilenameVectorProperties, "UTF-8"));
            outputFileVectorProperties.setIndent("    ");
        } catch (IOException e) {
            throw new UncheckedIOException("Unable to create output file", e);
        }
        int currSample = 0;
        int currExample = 0;
        int[] featuresDimensions = null;
        int[] labelsDimensions = null;
        ProgressLogger pg = new ProgressLogger(LOG);
        pg.expectedUpdates = baseInformationIterator.numExamples();
        pg.displayLocalSpeed = true;
        pg.itemsName = "records";
        pg.start();
        while (baseInformationIterator.hasNext()) {
            DataSet ds = baseInformationIterator.next();
            INDArray dsFeatures = ds.getFeatures();
            INDArray dsLabels = ds.getLabels();
            if (featuresDimensions == null) {
                featuresDimensions = dsFeatures.shape();
            }
            if (labelsDimensions == null) {
                labelsDimensions = dsLabels.shape();
            }
            writeIds(outputFileVector, currSample, currExample, 0);
            writeVectorFromINDArray(outputFileVector, dsFeatures);
            writeIds(outputFileVector, currSample, currExample, 1);
            writeVectorFromINDArray(outputFileVector, dsLabels);
            currExample++;
            pg.lightUpdate();
        }
        pg.stop();

        VectorProperties.VectorPropertiesSample[] samplesProps = new VectorProperties.VectorPropertiesSample[1];
        samplesProps[0] = new VectorProperties.VectorPropertiesSample(FilenameUtils.getBaseName(args().trainingSet),
                "test");
        VectorProperties.VectorPropertiesVector[] vectorsProps = new VectorProperties.VectorPropertiesVector[2];
        vectorsProps[0] = new VectorProperties.VectorPropertiesVector("features", "32_bit_float",
                featuresDimensions);
        vectorsProps[1] = new VectorProperties.VectorPropertiesVector("labels", "32_bit_float",
                labelsDimensions);
        VectorProperties vectorProperties = new VectorProperties(0, 1, samplesProps,
                vectorsProps);
        Gson gson = new Gson();
        gson.toJson(vectorProperties, VectorProperties.class, outputFileVectorProperties);
    }

    private void writeIds(PrintWriter writer, int sampleId, int exampleId, int vectorId) {
        writer.append(Integer.toString(sampleId))
                .append(" ")
                .append(Integer.toString(exampleId))
                .append(" ")
                .append(Integer.toString(vectorId));
    }

    private void writeVectorFromINDArray(PrintWriter writer, INDArray array) {
        NdIndexIterator ndIndexIterator = new NdIndexIterator('c', array.shape());
        while (ndIndexIterator.hasNext()) {
            writer.append(" ")
                .append(Float.toString(array.getFloat(ndIndexIterator.next())));
        }
        writer.append("\n");
    }

    static class VectorProperties {
        private int majorVersion;
        private int minorVersion;
        private VectorPropertiesSample[] samples;
        private VectorPropertiesVector[] vectors;

        VectorProperties(int majorVersion, int minorVersion, VectorPropertiesSample[] samples,
                                VectorPropertiesVector[] vectors) {
            this.majorVersion = majorVersion;
            this.minorVersion = minorVersion;
            this.samples = samples;
            this.vectors = vectors;
        }

        public int getMajorVersion() {
            return majorVersion;
        }

        public int getMinorVersion() {
            return minorVersion;
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
