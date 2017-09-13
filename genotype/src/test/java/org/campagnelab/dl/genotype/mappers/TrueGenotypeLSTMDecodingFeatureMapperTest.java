package org.campagnelab.dl.genotype.mappers;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.framework.mappers.MappedDimensions;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.Properties;

import static org.junit.Assert.assertEquals;

/**
 * Created by joshuacohen on 2/7/17.
 */
public class TrueGenotypeLSTMDecodingFeatureMapperTest {

    @Test
    public void testRecord1NotPredicting() throws Exception {
        String expectedInput = "[[[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask = "[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00]";
        TrueGenotypeLSTMDecodingFeatureMapper mapper = new TrueGenotypeLSTMDecodingFeatureMapper();
        Properties mapperProperties = new Properties();
        mapperProperties.setProperty("isPredicting", "false");
        mapper.configure(mapperProperties);
        MappedDimensions dim = mapper.dimensions();
        INDArray inputs = Nd4j.zeros(1, dim.numElements(1), dim.numElements(2));
        INDArray mask = Nd4j.zeros(1, dim.numElements(2));
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(TrueGenotypeLSTMLabelMapperTest.records[0], builder);
        final BaseInformationRecords.BaseInformation recordObject = builder.build();
        mapper.prepareToNormalize(recordObject, 0);
        mapper.mapFeatures(recordObject, inputs, 0);
        mapper.maskFeatures(recordObject, mask, 0);
        assertEquals(expectedInput, inputs.toString());
        assertEquals(expectedMask, mask.toString());
    }

    @Test
    public void testRecord2NotPredicting() throws Exception {
        String expectedLabel = "[[[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  1.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.00],  \n" +
                "  [1.00,  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask = "[1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00]";
        TrueGenotypeLSTMDecodingFeatureMapper mapper = new TrueGenotypeLSTMDecodingFeatureMapper();
        Properties mapperProperties = new Properties();
        mapperProperties.setProperty("isPredicting", "false");
        mapper.configure(mapperProperties);
        MappedDimensions dim = mapper.dimensions();
        INDArray inputs = Nd4j.zeros(1, dim.numElements(1), dim.numElements(2));
        INDArray mask = Nd4j.zeros(1, dim.numElements(2));
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(TrueGenotypeLSTMLabelMapperTest.records[1], builder);
        final BaseInformationRecords.BaseInformation recordObject = builder.build();
        mapper.prepareToNormalize(recordObject, 0);
        mapper.mapFeatures(recordObject, inputs, 0);
        mapper.maskFeatures(recordObject, mask, 0);
        assertEquals(expectedLabel, inputs.toString());
        assertEquals(expectedMask, mask.toString());
    }

    @Test
    public void testRecordsPredicting() throws Exception {
        String expectedInput = "[[[0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00],  \n" +
                "  [0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]]]";
        String expectedMask = "[1.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00]";
        TrueGenotypeLSTMDecodingFeatureMapper mapper = new TrueGenotypeLSTMDecodingFeatureMapper();
        Properties mapperProperties = new Properties();
        mapperProperties.setProperty("isPredicting", "true");
        mapper.configure(mapperProperties);
        MappedDimensions dim = mapper.dimensions();
        for (String record : TrueGenotypeLSTMLabelMapperTest.records) {
            INDArray inputs = Nd4j.zeros(1, dim.numElements(1), dim.numElements(2));
            INDArray mask = Nd4j.zeros(1, dim.numElements(2));
            final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
            TextFormat.getParser().merge(record, builder);
            final BaseInformationRecords.BaseInformation recordObject = builder.build();
            mapper.prepareToNormalize(recordObject, 0);
            mapper.mapFeatures(recordObject, inputs, 0);
            mapper.maskFeatures(recordObject, mask, 0);
            assertEquals(expectedInput, inputs.toString());
            assertEquals(expectedMask, mask.toString());
        }
    }
}