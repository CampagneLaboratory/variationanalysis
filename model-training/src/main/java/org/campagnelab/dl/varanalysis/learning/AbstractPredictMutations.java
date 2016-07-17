package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.calibrate.CalibratingModel;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.PrintWriter;
import java.util.List;
import java.util.SortedSet;

/**
 * Created by fac2003 on 6/10/16.
 */
public abstract class AbstractPredictMutations {

    final String header = "mutatedLabel\tProbabilityMut\tProbabilityUnmut\tcorrectness\tfrequency\tmutatedBase\trefIdx\tposition\treferenceBase\tsample1Counts\tsample2Counts\tsample1Scores\tsample2Scores\tsumCounts\tformatted1\tformatted2";

    protected void writeHeader(PrintWriter results) {
        results.append(header);
        if (cmodel!=null) {
            results.append("\tcalibratedP");
        }
        results.append("\n");
    }

    protected void writeRecordResult(MultiLayerNetwork model, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record, AreaUnderTheROCCurve aucLossCalculator) {
        writeRecordResult(model, null, results, featureMapper, pgReadWrite, record, aucLossCalculator, null, null);
    }

    CalibratingModel cmodel;

    protected void writeRecordResult(MultiLayerNetwork model, MultiLayerNetwork calibrationModel, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record, AreaUnderTheROCCurve aucLossCalculator, SortedSet<Float> plantedMutSet, SortedSet<Float> allRecSet) {

        INDArray testFeatures = Nd4j.zeros(1, featureMapper.numberOfFeatures());
        featureMapper.mapFeatures(record, testFeatures, 0);
        INDArray testPredicted = model.output(testFeatures, false);
        String features = featuresToString(record);
        //boolean
        boolean mutated = record.getMutated();
        ProtoPredictor predictor = new ProtoPredictor(model, featureMapper);
        ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
        String formatted0 = record.getSamples(0).getFormattedCounts().replaceAll("\n", "");
        String formatted1 = record.getSamples(1).getFormattedCounts().replaceAll("\n", "");
        String correctness = (prediction.clas == mutated) ? "right" : "wrong";
        if (aucLossCalculator != null) {
            aucLossCalculator.observe(prediction.posProb, mutated ? 1 : -1);
        }


        results.append(String.format("%s\t%f\t%f\t%s\t%s\t%s\t%s",
                (mutated ? "1" : "0"),
                prediction.posProb, prediction.negProb,
                correctness, features,
                formatted0, formatted1
        ));
        if (cmodel != null) {
            results.append(String.format("\t%f", cmodel.estimateCalibratedP(testFeatures)));
        }

        results.append("\n");
        pgReadWrite.update();

        //update tree sets
        if (plantedMutSet != null)

        {
            allRecSet.add(prediction.posProb);
            if (mutated) {
                plantedMutSet.add(prediction.posProb);
            }
        }

    }

    protected abstract String featuresToString(BaseInformationRecords.BaseInformation record);
}
