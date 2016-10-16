package org.campagnelab.dl.varanalysis.learning;

import org.apache.commons.collections.FastArrayList;
import org.apache.commons.compress.utils.IOUtils;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.campagnelab.dl.varanalysis.storage.RecordReader;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;

import java.io.IOException;
import java.util.List;

/**
 * Helper class to estimate performance on a test set, or part of it.
 * Created by fac2003 on 7/15/16.
 */
public class MeasurePerformance {
    private int aucClipMaxObservations;
    private int scoreN = Integer.MAX_VALUE;

    public MeasurePerformance(int scoreN, String datafilePath) throws IOException {
        this.scoreN = scoreN;
        this.aucClipMaxObservations = Integer.MAX_VALUE;
        readValidation(datafilePath);
    }

    public MeasurePerformance(int scoreN, int aucClipMaxObservations, String datafilePath) throws IOException {
        this(scoreN, datafilePath);
        this.aucClipMaxObservations = aucClipMaxObservations;
    }

    List<BaseInformationRecords.BaseInformation> records = new FastArrayList();

    private void readValidation(String datafilePath) throws IOException {
        int index = 0;
        RecordReader reader = new RecordReader(datafilePath);
        try {
            for (BaseInformationRecords.BaseInformation record : reader) {
                records.add(record);
                if (index > scoreN) break;
            }
        } finally {
            IOUtils.closeQuietly(reader);
        }
    }


    public double estimateAUC(FeatureMapper featureMapper, MultiLayerNetwork model) throws IOException {


        //may need to adjust batch size and write outputs piecewise if test sets are very large
        //BaseInformationIterator baseIter = new BaseInformationIterator(testsetPath, Integer.MAX_VALUE, new FeatureMapperV2(), new SimpleFeatureCalculator());

        //DataSet ds = baseIter.next();
//set up logger

        AreaUnderTheROCCurve aucLossCalculator = new AreaUnderTheROCCurve(aucClipMaxObservations);
        int index = 0;
        for (BaseInformationRecords.BaseInformation record : records) {
            boolean mutated = record.getMutated();
            ProtoPredictor predictor = new ProtoPredictor(model, featureMapper);
            ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
            aucLossCalculator.observe(prediction.posProb, mutated ? 1 : -1);

            index++;
            if (index > scoreN) break;
        }
        return aucLossCalculator.evaluateStatistic();


    }


}
