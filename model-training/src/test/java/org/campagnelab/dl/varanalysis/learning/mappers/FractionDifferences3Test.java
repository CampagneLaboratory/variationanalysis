package org.campagnelab.dl.varanalysis.learning.mappers;

import com.google.protobuf.TextFormat;
import org.apache.uima.util.FileUtils;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.FractionDifferences3;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.FileReader;

import static org.junit.Assert.assertNotEquals;

/**
 * Created by fac2003 on 6/11/16.
 */
public class FractionDifferences3Test {


    @Test
    public void mapFeatures() throws Exception {
        int index = 0;

        for (String[] recordFilenamePair : recordsFilenames) {

            String recordFilename_plain = recordFilenamePair[0];
            String recordFilename_mutated = recordFilenamePair[1];
            String record = FileUtils.reader2String(new FileReader(recordFilename_plain));
            String recordMutated = FileUtils.reader2String(new FileReader(recordFilename_mutated));


            final BaseInformationRecords.BaseInformation.Builder builderPlain = getBuilder(record);
            final BaseInformationRecords.BaseInformation.Builder builderMutated = getBuilder(recordMutated);

            INDArray inputs_plain = getFeatures(builderPlain);
            INDArray inputs_mutated = getFeatures(builderMutated);
            System.out.println("plain:    "+inputs_plain.toString());
            System.out.println("mutated:  "+inputs_mutated.toString());
            assertNotEquals("features for plain and mutated must differ", inputs_plain.toString(), inputs_mutated.toString());

            index++;
        }
    }

    private INDArray getFeatures(BaseInformationRecords.BaseInformation.Builder builderPlain) {
        FeatureMapper calculator = new FractionDifferences3();
        INDArray inputs = Nd4j.zeros(1, calculator.numberOfFeatures());
        calculator.prepareToNormalize(builderPlain.build(), 0);
        calculator.mapFeatures(builderPlain.build(), inputs, 0);
        return inputs;
    }

    private BaseInformationRecords.BaseInformation.Builder getBuilder(String record) throws TextFormat.ParseException {
        final BaseInformationRecords.BaseInformation.Builder builder = BaseInformationRecords.BaseInformation.newBuilder();
        TextFormat.getParser().merge(record, builder);
        return builder;
    }

    String[][] recordsFilenames = {{"test-data/focus-on-errors/11_9905397.txt", "test-data/focus-on-errors/11_9905397_mutated.txt"}};
    String[] expectedFeatures = {""};
}