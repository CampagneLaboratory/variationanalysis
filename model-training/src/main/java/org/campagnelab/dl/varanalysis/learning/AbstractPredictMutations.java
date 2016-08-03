package org.campagnelab.dl.varanalysis.learning;

import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.model.utils.CalcCalibrator;
import org.campagnelab.dl.model.utils.ProtoPredictor;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.learning.calibrate.CalibratingModel;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.stats.AreaUnderTheROCCurve;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;

/**
 * Created by fac2003 on 6/10/16.
 */
public abstract class AbstractPredictMutations {

    String header;
    final String s3Scores = "\tsample3Scores";
    final String s3Counts = "\tsampleS3Counts";
    final String formatted3 = "\tformatted3";


    protected void writeHeader(PrintWriter results) {

        header = "mutatedLabel\tProbabilityMut\tProbabilityUnmut\tcorrectness\tfrequency" +
                "\tmutatedBase\trefIdx\tposition\treferenceBase" +
                "\tsample1Counts\tsample2Counts" + s3Counts +
                "\tsample1Scores\tsample2Scores" + s3Scores +
                "\tsumCounts\tformatted1\tformatted2" + formatted3;
        results.append(header);
        if (cmodel!=null) {
            results.append("\tcalibratedP");
        }
        results.append("\n");
    }

    protected void writeRecordResult(MultiLayerNetwork model, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record, AreaUnderTheROCCurve aucLossCalculator, boolean isTrio) {
        writeRecordResult(model, null, results, featureMapper, pgReadWrite, record, aucLossCalculator, null, isTrio);
    }
    protected void writeRecordResult(MultiLayerNetwork model, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record, AreaUnderTheROCCurve aucLossCalculator) {
        writeRecordResult(model, null, results, featureMapper, pgReadWrite, record, aucLossCalculator, null, false);
    }

    CalibratingModel cmodel;

    protected void writeRecordResult(MultiLayerNetwork model, MultiLayerNetwork calibrationModel, PrintWriter results, FeatureMapper featureMapper, ProgressLogger pgReadWrite, BaseInformationRecords.BaseInformation record, AreaUnderTheROCCurve aucLossCalculator, CalcCalibrator calc, boolean isTrio) {
        INDArray testFeatures = Nd4j.zeros(1, featureMapper.numberOfFeatures());
        featureMapper.mapFeatures(record, testFeatures, 0);
        INDArray testPredicted = model.output(testFeatures, false);
        String features = featuresToString(record);
        //boolean
        boolean mutated = record.getMutated();
        ProtoPredictor predictor = new ProtoPredictor(model, featureMapper);
        ProtoPredictor.Prediction prediction = predictor.mutPrediction(record);
        String formatted0 = genFormattedString(record.getSamples(0));
        String formatted1 = genFormattedString(record.getSamples(1));
        String formatted2 = isTrio?"\t"+genFormattedString(record.getSamples(2)):"";
        String correctness = (prediction.clas == mutated) ? "right" : "wrong";
        if (aucLossCalculator != null) {
            aucLossCalculator.observe(prediction.posProb, mutated ? 1 : -1);
        }


        results.append(String.format("%s\t%f\t%f\t%s\t%s\t%s\t%s%s",
                (mutated ? "1" : "0"),
                prediction.posProb, prediction.negProb,
                correctness, features,
                formatted0, formatted1, formatted3
        ));
        if (cmodel != null) {
            results.append(String.format("\t%f", cmodel.estimateCalibratedP(testFeatures)));
        }

        results.append("\n");
        pgReadWrite.update();

        //update tree sets
        calc.observe(prediction.posProb,mutated);

    }

    private String genFormattedString(BaseInformationRecords.SampleInfo sample) {
        int a = sample.getCounts(0).getGenotypeCountReverseStrand() + sample.getCounts(0).getGenotypeCountForwardStrand();
        int t = sample.getCounts(1).getGenotypeCountReverseStrand() + sample.getCounts(1).getGenotypeCountForwardStrand();
        int c = sample.getCounts(2).getGenotypeCountReverseStrand() + sample.getCounts(2).getGenotypeCountForwardStrand();
        int g = sample.getCounts(3).getGenotypeCountReverseStrand() + sample.getCounts(3).getGenotypeCountForwardStrand();
        int n = sample.getCounts(4).getGenotypeCountReverseStrand() + sample.getCounts(4).getGenotypeCountForwardStrand();
        String fb = sample.getFormattedCounts().split(" ")[8];
        int numIndels = sample.getCountsCount() - 5;
        int[] indels = new int[numIndels];
        for (int i = 5; i < numIndels + 5 ; i++) {
            indels[i-5] = sample.getCounts(i).getGenotypeCountForwardStrand() + sample.getCounts(i).getGenotypeCountReverseStrand();
        }
        return String.format("counts A=%d T=%d C=%d G=%d N=%d %s indels:%s",a,t,c,g,n,fb, Arrays.toString(indels));
    }

    protected abstract String featuresToString(BaseInformationRecords.BaseInformation record);
}
