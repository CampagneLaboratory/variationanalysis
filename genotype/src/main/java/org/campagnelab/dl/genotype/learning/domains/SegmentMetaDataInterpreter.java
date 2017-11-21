package org.campagnelab.dl.genotype.learning.domains;

import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.predictions.SegmentMetaData;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * Created by fac2003 on 11/1/17.
 */
public class SegmentMetaDataInterpreter implements PredictionInterpreter<SegmentInformationRecords.SegmentInformation,
        SegmentMetaData> {

    @Override
    public SegmentMetaData interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        SegmentMetaData metaData = new SegmentMetaData();
        metaData.index = exampleIndex;
        //final int numLabels = output.size(1);
        final int sequenceLength = output.size(2);
        INDArray trueLabelRow = trueLabels.getRow(exampleIndex);

        for (int baseIndex = 0; baseIndex < sequenceLength; baseIndex++) {


            double encodedMetadata = trueLabelRow.getDouble(baseIndex);
            if (encodedMetadata > 0.1 && encodedMetadata < 0.9) {
                metaData.isSnp.set(baseIndex);
            } else if (encodedMetadata > 1.1) {
                metaData.isIndel.set(baseIndex);
            } else {
                metaData.isSnp.clear(baseIndex);
                metaData.isIndel.clear(baseIndex);
            }
        }
        return metaData;
    }

    @Override
    public SegmentMetaData interpret(SegmentInformationRecords.SegmentInformation record, INDArray output) {
        SegmentMetaData metaData = new SegmentMetaData();
        SegmentInformationRecords.Sample sample = record.getSample(0);
        for (int baseIndex = 0; baseIndex < record.getLength(); baseIndex++) {
            if (sample.getBase(baseIndex).getIsVariant()) {
                metaData.isSnp.set(baseIndex);
            }
            if (sample.getBase(baseIndex).getHasTrueIndel()) {
                metaData.isIndel.set(baseIndex);
            }
        }
        return metaData;
    }
}
