package org.campagnelab.dl.varanalysis.learning.mappers;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.FeatureMapper;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.PrintWriter;

/**
 * Created by fac2003 on 6/10/16.
 */
public abstract class AbstractPredictMutations {

    final String header = "mutatedLabel\tProbabilityMut\tProbabilityUnmut\tcorrectness\tfrequency\tmutatedBase\trefIdx\tposition\treferenceBase\tsample1Counts\tsample2Counts\tsample1Scores\tsample2Scores\tsumCounts\tformatted1\tformatted2\n";

    protected void writeHeader(PrintWriter results) {
        results.append(header);
    }

    protected void writeRecordResult(MultiLayerNetwork model, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record) {
        INDArray testFeatures = Nd4j.zeros(1, featureMapper.numberOfFeatures());
        featureMapper.mapFeatures(record,testFeatures,0);
        INDArray testPredicted = model.output(testFeatures,false);
        String features = featuresToString(record);
        //boolean
        boolean mutated = record.getMutated();
        ProtoPredictor predictor = new ProtoPredictor(model,featureMapper);
        ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
        String formatted0 = record.getSamples(0).getFormattedCounts().replaceAll("\n","");
        String formatted1 = record.getSamples(1).getFormattedCounts().replaceAll("\n","");
        String correctness = (prediction.clas == mutated) ? "right" : "wrong";
        results.append((mutated?"1":"0") + "\t" + Float.toString(prediction.posProb) + "\t" + Float.toString(prediction.negProb) +  "\t" + correctness + "\t" + features + "\t" + formatted0 + "\t" + formatted1 + "\n");
        pgReadWrite.update();
    }

    protected abstract String featuresToString(BaseInformationRecords.BaseInformation record);
}
