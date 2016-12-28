package org.campagnelab.dl.framework.mappers;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Maps a one dimensional feature of numFeatures ones or zeros, with configuration over
 * whether the feature should be masked or not
 *
 * Created by joshuacohen on 12/14/16.
 */
public class NAryFeatureMapper<RecordType> implements FeatureMapper<RecordType> {
    static private Logger LOG = LoggerFactory.getLogger(NAryFeatureMapper.class);

    private int numFeatures;
    private boolean isZero;
    private boolean hasMask;
    private boolean isMasked;
    private MappedDimensions dim;

    private int[] mapperIndices = new int[]{0, 0};

    /**
     * @param numFeatures Number of individual features the feature mapper maps
     * @param isZero If true, feature filled with zeros; if false, filled with ones
     * @param hasMask Whether or not the mapper has a mask
     * @param isMasked If hasMask is true, specifies if mask should be 1 (i.e., isMasked true) or 0
     */
    public NAryFeatureMapper(int numFeatures, boolean isZero, boolean hasMask, boolean isMasked) {
        this.numFeatures = numFeatures;
        this.isZero = isZero;
        if (numFeatures > 1 && !isZero) {
            LOG.warn("Mapping multiple 1s may break one-hot encoding");
        }
        this.hasMask = hasMask;
        this.isMasked = isMasked;
        dim = new MappedDimensions(numFeatures);
    }

    @Override
    public int numberOfFeatures() {
        return numFeatures;
    }

    @Override
    public MappedDimensions dimensions() {
        return dim;
    }

    @Override
    public void prepareToNormalize(RecordType record, int indexOfRecord) {

    }

    @Override
    public void mapFeatures(RecordType record, INDArray inputs, int indexOfRecord) {
        mapperIndices[0] = indexOfRecord;
        for (int i = 0; i < numFeatures; i++) {
            mapperIndices[1] = i;
            inputs.putScalar(mapperIndices, produceFeature(record, i));
        }
    }

    @Override
    public boolean hasMask() {
        return hasMask;
    }

    @Override
    public void maskFeatures(RecordType record, INDArray mask, int indexOfRecord) {
        mask.putScalar(indexOfRecord, isMasked(record, -1) ? 1F : 0F);
    }

    @Override
    public boolean isMasked(RecordType record, int featureIndex) {
        return isMasked;
    }

    @Override
    public float produceFeature(RecordType record, int featureIndex) {
        return isZero ? 0F : 1F;
    }
}
