package org.campagnelab.dl.varanalysis.learning;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.models.ModelLoader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;

import java.io.IOException;
import java.util.Properties;

import static org.junit.Assert.assertTrue;

/**
 * Created by fac2003 on 7/15/16.
 */
public class MeasurePerformanceTest {
    //@Test
    public void testEstimateAUC() throws IOException {
        MeasurePerformance perf=new MeasurePerformance(1000,"/data/mutated-MHFC-13-CTL_B_NK.parquet");
     ModelLoader loader=new ModelLoader("sample_data/1468533491636");
        FeatureMapper mapper=loader.loadFeatureMapper(new Properties());
        MultiLayerNetwork model = loader.loadModel("best");
        double auc=perf.estimateAUC(mapper, model );
        assertTrue(auc>0.5);
    }
}