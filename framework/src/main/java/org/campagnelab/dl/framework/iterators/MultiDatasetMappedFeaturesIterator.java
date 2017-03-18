package org.campagnelab.dl.framework.iterators;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Properties;

/**
 * An Iterator over mapped features (.cf/.cfp files).
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/2/16.
 */
public class MultiDatasetMappedFeaturesIterator implements MultiDataSetIterator {
    static private Logger LOG = LoggerFactory.getLogger(MappedFeaturesIterator.class);

    private final int numExamples;
    private final int miniBatchSize;
    private final FastBufferedInputStream inputStream;
    private final int cacheN;
    private int index;
    private MultiDataSetPreProcessor preProcessor;
    private int numDevices;

    public MultiDatasetMappedFeaturesIterator(String basename) {
        this(basename, Integer.MAX_VALUE);
    }

    public MultiDatasetMappedFeaturesIterator(String basename, int cacheN) {
        try {
            Properties cfProperties = new Properties();
            cfProperties.load(new FileReader(basename + ".cfp"));

            miniBatchSize = Integer.parseInt(cfProperties.getProperty("miniBatchSize", "0"));
            numExamples = Integer.parseInt(cfProperties.getProperty("numRecords", "0"));
            inputStream = new FastBufferedInputStream(new FileInputStream(basename + ".cf"));
            this.cacheN = cacheN;
        } catch (Exception e) {
            throw new RuntimeException("Unable to create MappedFeaturesIterator ", e);
        }
    }

    @Override
    public MultiDataSet next(int miniBatchSize) {
        if (miniBatchSize != this.miniBatchSize) {
            throw new IllegalArgumentException("numExamples must match the cached minibatchSize: " + miniBatchSize);
        }
        return next();
    }

    @Override
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor) {
        this.preProcessor = preProcessor;
    }

    @Override
    public boolean resetSupported() {
        return true;
    }


    public boolean asyncSupported() {
        return true;
    }

    @Override
    public void reset() {
        index = 0;
        try {
            inputStream.position(0);
        } catch (IOException e) {
            LOG.error("Unable to reset iterator to position 0");
        }
        numDevices = Nd4j.getAffinityManager().getNumberOfDevices();

    }


    @Override
    public boolean hasNext() {

        return index < Math.min(numExamples, cacheN);
    }

    byte[] length = new byte[4];
    ByteArrayList content = new ByteArrayList();
    long i = 0;

    @Override
    public MultiDataSet next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        try {
            inputStream.read(length, 0, 4);
        } catch (IOException e) {
            LOG.error("Unable to read length from stream.", e);
        }
        int l = (length[0] << 8 * 3 & 0xFF000000) |
                (length[1] << 8 * 2 & 0x00FF0000) |
                (length[2] << 8 & 0x0000FF00) |
                (length[3] & 0x000000FF);
        // resize the buffer if needed:
        content.size(l);
        final byte[] elements = content.elements();
        try {
            inputStream.read(elements);
        } catch (IOException e) {
            LOG.error("Unable to read content from stream at index " + index, e);
        }
        MultiDataSet ds = new org.nd4j.linalg.dataset.MultiDataSet();
        if (numDevices >= 2) {
            // try to move the MDS to a random device, so that we won't run out of memory on the first one
            // as would happen if all MDS went there (they were built on a cpu)

            moveToDecide(ds, (int) (i++ % numDevices));
        }

        try (ByteArrayInputStream from = new ByteArrayInputStream(elements)) {
            ds.load(from);
        } catch (IOException e) {
            LOG.error("Unable to load dataset at index " + index, e);
        }
        if (preProcessor != null) {
            preProcessor.preProcess(ds);
        }
        index += ds.getFeatures(0).size(0);
        return ds;

    }

    private void moveToDecide(MultiDataSet multiDataSet, int deviceIndex) {
        for (INDArray array : multiDataSet.getFeatures()) {
            Nd4j.getAffinityManager().replicateToDevice(deviceIndex, array);
        }
        for (INDArray array : multiDataSet.getLabels()) {
            Nd4j.getAffinityManager().replicateToDevice(deviceIndex, array);
        }
        for (INDArray array : multiDataSet.getLabelsMaskArrays()) {
            Nd4j.getAffinityManager().replicateToDevice(deviceIndex, array);
        }
        for (INDArray array : multiDataSet.getLabelsMaskArrays()) {
            Nd4j.getAffinityManager().replicateToDevice(deviceIndex, array);
        }
    }
}

