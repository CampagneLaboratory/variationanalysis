package org.campagnelab.dl.framework.mixup;

import cern.jet.random.Beta;
import cern.jet.random.engine.RandomEngine;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;

/**
 * Implements mixup as a preprocessor on a multi-dataset.
 * Created by fac2003 on 11/1/17.
 */
public class MixupMultiDataSetPreProcessor implements MultiDataSetPreProcessor {
    private final double alpha;
    XoRoShiRo128PlusRandom random;

    Beta beta;

    public MixupMultiDataSetPreProcessor(long seed, double alpha) {
        random = new XoRoShiRo128PlusRandom(seed);
        assert  alpha>0:"alpha must be strictly positive.";
        beta = new Beta(alpha, alpha, new RandomEngine() {
            @Override
            public int nextInt() {
                return random.nextInt();
            }
        });
        this.alpha = alpha;
    }

    @Override
    public void preProcess(MultiDataSet multiDataSet) {
        double alm = beta.nextDouble();
        INDArray[] features = multiDataSet.getFeatures();
        INDArray[] labels = multiDataSet.getLabels();
        INDArray[] featureMasks = multiDataSet.getFeaturesMaskArrays();
        INDArray[] labelMasks = multiDataSet.getLabelsMaskArrays();


        int minibatchSize = features[0].size(0);
        INDArray[] tmpBuffer = new INDArray[minibatchSize];
        int randomIndex1[] = new int[minibatchSize];
        int randomIndex2[] = new int[minibatchSize];
        // determine how this minibatch will be mixuped:
        for (int exampleIndex = 0; exampleIndex < minibatchSize; exampleIndex++) {

            tmpBuffer[exampleIndex] = null;// Nd4j.zeros(features[inputIndex].shape());
            // draw two indices randomly from this minibatch, keep track of these indices
            // to mix labels and mask appropriately.
            randomIndex1[exampleIndex] = random.nextInt(minibatchSize);
            randomIndex2[exampleIndex] = random.nextInt(minibatchSize);
        }
        for (INDArray feature : features) shuffle(minibatchSize, alm, feature, tmpBuffer, randomIndex1, randomIndex2);

        for (INDArray label : labels) shuffle(minibatchSize, alm, label, tmpBuffer, randomIndex1, randomIndex2);

        for (INDArray featureMask : featureMasks)
            keepLongestMask(minibatchSize, featureMask, tmpBuffer, randomIndex1, randomIndex2);

        for (INDArray labelMask : labelMasks)
            keepLongestMask(minibatchSize, labelMask, tmpBuffer, randomIndex1, randomIndex2);
    }

    private void keepLongestMask(int minibatchSize,  INDArray mask, INDArray[] tmpBuffer, int[] randomIndex1, int[] randomIndex2) {
     if (mask==null ) return;
        // Find the longest mask and keep it as mixup ask:
        for (int exampleIndex = 0; exampleIndex < minibatchSize; exampleIndex++) {
            int random1 = randomIndex1[exampleIndex];
            int random2 = randomIndex2[exampleIndex];
            final INDArray mask1 = mask.getRow(random1);
            final INDArray mask2 = mask.getRow(random2);
            if (mask1.sub(mask2).sumNumber().doubleValue() < 0) {
                // mask2 has more 1s than mask1, use mask2:
                tmpBuffer[exampleIndex] = mask2;
            } else {
                tmpBuffer[exampleIndex] = mask1;
            }

        }
        for (int exampleIndex = 0; exampleIndex < minibatchSize; exampleIndex++) {

            // assign tmpBuffer[inputIndex] back into the minibatch:
            mask.putRow(exampleIndex, tmpBuffer[exampleIndex]);
        }
    }

    private void shuffle(int minibatchSize, double alm, INDArray features, INDArray[] tmpBuffer, int[] randomIndex1, int[] randomIndex2) {
        for (int exampleIndex = 0; exampleIndex < minibatchSize; exampleIndex++) {
            int random1 = randomIndex1[exampleIndex];
            int random2 = randomIndex2[exampleIndex];
            final INDArray example1 = features.getRow(random1);
            final INDArray example2 = features.getRow(random2);
            // new example is linear combination of example 1 and example2:
            tmpBuffer[exampleIndex] = example1.mul(alm).addi(example2.mul(1.0 - alm));
        }
        for (int exampleIndex = 0; exampleIndex < minibatchSize; exampleIndex++) {

            // assign tmpBuffer[inputIndex] back into the minibatch:
            features.putRow(exampleIndex, tmpBuffer[exampleIndex]);
        }

    }


}
