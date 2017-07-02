package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.genotype.learning.domains.predictions.SingleGenotypePrediction;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.List;

/**
 * Decoder for SoftmaxGenotypePrediction.
 * Created by fac2003 on 7/2/17.
 */
public class AggregatedSoftmaxGenotypePrediction extends GenotypePrediction {
    ObjectArraySet<String> predictedAlleles;
    ObjectArraySet<String> trueAlleles;

    public AggregatedSoftmaxGenotypePrediction(BaseInformationRecords.BaseInformation record,
                                               List<Prediction> individualOutputPredictions) {
        predictedAlleles = new ObjectArraySet<>();
        trueAlleles = new ObjectArraySet<>();
        set(record, (SoftmaxGenotypePrediction) individualOutputPredictions.get(0),
                (MetadataPrediction) individualOutputPredictions.get(1));
    }

    private void set(BaseInformationRecords.BaseInformation record,
                     SoftmaxGenotypePrediction softmaxGenotype, MetadataPrediction metaData) {
        int numBits = softmaxGenotype.numBits;
        for (int i = 0; i < numBits; i++) {
            // we decode the sorted allele index:
            boolean alleleIsCalled = ((softmaxGenotype.predictedGenotypeIndex >> i) & 0x1) != 0;
            boolean alleleIsInTrueGenotype = ((softmaxGenotype.trueGenotypeIndex >> i) & 0x1) != 0;
            if (alleleIsCalled || alleleIsInTrueGenotype) {
                SingleGenotypePrediction pred = new SingleGenotypePrediction();
                pred.probabilityIsCalled = softmaxGenotype.probability;
                pred.sortedCountIndex = i;
                pred.trueIsCalled = alleleIsInTrueGenotype;

                String allele = pred.predictedSingleGenotype(record, metaData);
                if (alleleIsCalled) predictedAlleles.add(allele);
                if (alleleIsInTrueGenotype) trueAlleles.add(allele);
            }

        }
        set(metaData);
        this.overallProbability = softmaxGenotype.probability;
        this.isVariantProbability = 0.001;
        this.trueGenotype= GenotypeHelper.fromAlleles(trueAlleles);
        this.predictedGenotype= GenotypeHelper.fromAlleles(predictedAlleles);

        int value=1;
        value+=1;

    }
}
