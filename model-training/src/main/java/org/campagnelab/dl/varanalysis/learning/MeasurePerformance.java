package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Helper class to estimate performance on a test set, or part of it.
 * Created by fac2003 on 7/15/16.
 */
public class MeasurePerformance {
    int scoreN = Integer.MAX_VALUE;

    public MeasurePerformance(int scoreN) {
        this.scoreN = scoreN;
    }

    public double estimateAUC(FeatureMapper featureMapper, MultiLayerNetwork model, String datafilePath) throws IOException {


        //may need to adjust batch size and write outputs piecewise if test sets are very large
        //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());
        RecordReader reader = new RecordReader(datafilePath);
        //DataSet ds = baseIter.next();
//set up logger
        try {
            AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve();
            int index = 0;
            for (BaseInformationRecords.BaseInformation record : reader) {
                boolean mutated = record.getMutated();
                ProtoPredictor predictor = new ProtoPredictor(model, featureMapper);
                ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
                aucLossCalculator.observe(prediction.posProb, mutated ? 1 : -1);

                index++;
                if (index > scoreN) break;
            }
            return aucLossCalculator.evaluateStatistic();
        } finally {
            IOUtils.closeQuietly(reader);

        }

    }


}
