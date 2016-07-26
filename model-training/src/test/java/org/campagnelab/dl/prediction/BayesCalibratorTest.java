package org.campagnelab.dl.prediction;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.model.utils.BayesCalibrator;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import java.io.File;
import java.io.StringWriter;

import static org.junit.Assert.*;

/**
 * Test the mutator on some specific examples.
 * Created by fac2003 on 5/27/16.
 */
public class BayesCalibratorTest {

    @Test
    public void bayesTest() throws Exception {
        String modelDir = "test-results/testModel";
        new File(modelDir).mkdirs();
        BayesCalibrator calc = new BayesCalibrator(modelDir,"best",false);
        calc.observe(0.1f,false);
        calc.observe(0.2f,false);
        calc.observe(0.3f,false);
        calc.observe(0.4f,true);
        calc.observe(0.5f,true);
        calc.observe(0.6f,false);
        calc.observe(0.7f,true);
        calc.observe(0.8f,true);
        calc.observe(0.9f,false);
        calc.observe(1.0f,true);

        calc.save();

        calc = new BayesCalibrator(modelDir,"best",true);


        float expected = (1f/5f)/((4f*1000000f)*(1f/5f));
        System.out.println("expected:"+Float.toString(expected));
        assertEquals(expected,calc.calibrateProb(0.85f),0.000001f);


    };

}