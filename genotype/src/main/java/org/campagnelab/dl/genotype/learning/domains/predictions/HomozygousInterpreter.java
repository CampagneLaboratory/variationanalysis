package org.campagnelab.dl.genotype.learning.domains.predictions;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Interpret model output to determine if a site is homozygous/het and what hte homosygous genotype is.
 * Created by rct66 on 11/12/16.
 */
public class HomozygousInterpreter implements PredictionInterpreter<BaseInformationRecords.BaseInformation, HomozygousPrediction> {

    double maxProbability;
    boolean isHomozygous;
    private int maxIndex;

    public HomozygousInterpreter() {
    }


    @Override
    public HomozygousPrediction interpret(INDArray trueLabels, INDArray[] outputs, int predictionIndex) {
        throw new RuntimeException("a wrong interpret method was called on the homozygous interpeter");
    }

    @Override
    public HomozygousPrediction interpret(BaseInformationRecords.BaseInformation record, INDArray output) {
        HomozygousPrediction pred = new HomozygousPrediction();
        pred.inspectRecord(record);
        pred.predictedHomozygousGenotype = getHomozygousPrediction(record, output);
        pred.probability = maxProbability;
        pred.isHomozygous =(maxIndex!=10); // max probability on the isHomozygous site.
        return pred;
    }

    public String getHomozygousPrediction(BaseInformationRecords.BaseInformation currentRecord, INDArray output) {
        StringBuffer genotype = new StringBuffer();
        maxProbability = -1;
        maxIndex = -1;
        int predictionIndex = 0;
        for (int i = 0; i < 11; i++) {
            if (maxProbability < output.getDouble(predictionIndex, i)) {
                maxIndex = i;
                maxProbability = output.getDouble(predictionIndex, i);
            }
        }
        if (maxIndex > 10 || maxIndex == -1) {
            return "";
        }
        try {
            return currentRecord.getSamples(0).getCounts(maxIndex).getToSequence();
        } catch (IndexOutOfBoundsException e) {
            // predicted, but not present in the input features?
            return ".";
        }
    }
}
