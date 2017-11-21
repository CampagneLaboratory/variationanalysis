package org.campagnelab.dl.genotype.predictions;

import com.google.common.base.Joiner;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.genotype.mappers.SingleBaseLabelMapperV1;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class SegmentPredictionInterpreter implements PredictionInterpreter<SegmentInformationRecords.SegmentInformation, SegmentGenotypePrediction> {
    Int2ObjectMap<String> indicesToGenotypesMap;

    SingleBaseLabelMapperV1 mapper = new SingleBaseLabelMapperV1(0);

    /**
     * Create a SegmentPredictionInterpreter to interpret the "genotype" output.
     *
     * @param segmentPropertiesFilename Path the .ssip file where the mapping between genotypes and integer is defined.
     */
    public SegmentPredictionInterpreter(String segmentPropertiesFilename) {
        this.indicesToGenotypesMap = new Int2ObjectOpenHashMap<>();
        Properties ssip = new Properties();
        try {
            ssip.load(new FileInputStream(new File(segmentPropertiesFilename)));
            String numString = ssip.getProperty("genotype.segment.label.numOfEntries");
            assert numString != null : "property genotype.segment.label.numOfEntries must exist in .ssip";
            int numGenotypes = Integer.parseInt(numString);
            for (int i = 0; i < numGenotypes; i++) {
                String genotype = ssip.getProperty(String.format("genotype.segment.label.%d", i));
                indicesToGenotypesMap.put(i + 1, genotype);
            }
        } catch (IOException e) {
            throw new InternalError("Unable to load genotype mapping from ssip: " + segmentPropertiesFilename);
        }
    }

    @Override
    /**
     * Read genotypes and their probability from the output called "genotype".
     */
    public SegmentGenotypePrediction interpret(INDArray trueLabels, INDArray output, int exampleIndex) {
        SegmentGenotypePrediction prediction = new SegmentGenotypePrediction();
        prediction.index = exampleIndex;
        //final int numLabels = output.size(1);
        final int sequenceLength = output.size(1);
        prediction.probabilities = new float[sequenceLength];
        INDArray predictedRow = output;
        INDArray trueMaxIndices = null;
        if (trueLabels != null) {
            INDArray trueLabelRow = trueLabels;
            trueMaxIndices = Nd4j.argMax(trueLabelRow, 0);
        }
        INDArray predictedMaxIndices = Nd4j.argMax(predictedRow, 0);

        prediction.predictedGenotypes = new String[sequenceLength];
        prediction.trueGenotypes = new String[sequenceLength];
        for (int baseIndex = 0; baseIndex < sequenceLength; baseIndex++) {

            final int predictedMaxIndex = predictedMaxIndices.getInt(baseIndex);
            final int trueMaxIndex = (trueMaxIndices != null)?trueMaxIndices.getInt(baseIndex):0;
            prediction.probabilities[baseIndex] = predictedRow.getFloat(predictedMaxIndex);
            prediction.predictedGenotypes[baseIndex] = predictedMaxIndex == 0 ? null : indicesToGenotypesMap.get(predictedMaxIndex);
            prediction.trueGenotypes[baseIndex] = trueMaxIndex == 0 ? null : indicesToGenotypesMap.get(trueMaxIndex);
           // assume we know the length of the sequence (we do, since it is the number of features we collected.)
        }
        return prediction;
    }

    @Override
    public SegmentGenotypePrediction interpret(SegmentInformationRecords.SegmentInformation record, INDArray output) {
        SegmentGenotypePrediction prediction = interpret(null,output,0);
        prediction.length=record.getSample(0).getBaseCount();
        //loop on bases, and collect labels and populate the trueGenotypes in the prediction
       /* int numOfBases = calculateNumOfBases(record);
        prediction.probabilities = new float[numOfBases];
        prediction.trueGenotypes = new String[numOfBases];
        prediction.predictedGenotypes = new String[numOfBases];
        int baseIndex = 0;
        for (SegmentInformationRecords.Sample sample : record.getSampleList()) {
            for (SegmentInformationRecords.Base base : sample.getBaseList()) {
                INDArray predictedRow = output.getRow(baseIndex);
              

                prediction.probabilities[baseIndex] = predictedRow.getFloat(baseIndex);
                prediction.predictedGenotypes[baseIndex] = Joiner.on("/").join(base.getTrueLabelList());
                //prediction.trueGenotypes[baseIndex] = trueMaxIndex == 0 ? null : indicesToGenotypesMap.get(trueMaxIndex);
                if (prediction.trueGenotypes[baseIndex] != null) {
                    prediction.length += 1;
                }
                else {
                    break;
                }
                baseIndex++;
            }

        }   */

        return prediction;
    }

    private int calculateNumOfBases(SegmentInformationRecords.SegmentInformation record) {
        int b = 0;
        for (SegmentInformationRecords.Sample sample : record.getSampleList()) {
            for (SegmentInformationRecords.Base base : sample.getBaseList()) {
                  b++;
            }
        }
        return b;
    }
}
