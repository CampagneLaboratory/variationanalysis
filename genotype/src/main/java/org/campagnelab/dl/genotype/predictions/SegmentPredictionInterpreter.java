package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import org.campagnelab.dl.framework.domains.prediction.PredictionInterpreter;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class SegmentPredictionInterpreter implements PredictionInterpreter<SegmentInformationRecords.SegmentInformation, SegmentGenotypePrediction> {
    Int2ObjectMap<String> indicesToGenotypesMap;

    /**
     * Create a SegmentPredictionInterpreter to interpret the "genotype" output.
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
                String genotype=ssip.getProperty(String.format("genotype.segment.label.%d",i));
                indicesToGenotypesMap.put(i,genotype);
            }
        } catch (IOException e) {
            throw new InternalError("Unable to load genotype mapping from ssip: " + segmentPropertiesFilename);
        }
    }

    @Override
    /**
     * Read genotypes and their probability from the output called "genotype".
     */
    public SegmentGenotypePrediction interpret(INDArray trueLabels, INDArray output, int predictionIndex) {
        SegmentGenotypePrediction prediction = new SegmentGenotypePrediction();
        prediction.index = predictionIndex;
        //final int numLabels = output.size(1);
        final int sequenceLength=output.size(2);
        prediction.probabilities = new float[sequenceLength];
        INDArray predictedRow = output.getRow(predictionIndex);
        INDArray trueLabelRow = trueLabels.getRow(predictionIndex);
        INDArray predictedMaxIndices = Nd4j.argMax(predictedRow, 0);
        INDArray trueMaxIndices = Nd4j.argMax(trueLabelRow, 0);
        prediction.predictedGenotypes=new String[sequenceLength];
        prediction.trueGenotypes=new String[sequenceLength];
        for (int baseIndex = 0; baseIndex < sequenceLength; baseIndex++) {

            final int predictedMaxIndex = predictedMaxIndices.getInt(baseIndex);
            final int trueMaxIndex = trueMaxIndices.getInt(baseIndex);
            prediction.probabilities[baseIndex] = predictedRow.getFloat(predictedMaxIndex);
            prediction.predictedGenotypes[baseIndex] = indicesToGenotypesMap.get(predictedMaxIndex);
            prediction.trueGenotypes[baseIndex]=indicesToGenotypesMap.get(trueMaxIndex);
        }
        return prediction;
    }

    @Override
    public SegmentGenotypePrediction interpret(SegmentInformationRecords.SegmentInformation record, INDArray output) {
        throw new InternalError("Not currently implemented. Implement when PredictGS is needed.");
    }


}
