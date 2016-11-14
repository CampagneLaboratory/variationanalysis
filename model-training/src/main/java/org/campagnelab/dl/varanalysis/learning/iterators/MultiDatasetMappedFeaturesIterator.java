package org.campagnelab.dl.varanalysis.learning.iterators;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.lang.NotImplementedException;
import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.model.utils.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.domains.DomainDescriptor;
import org.campagnelab.dl.varanalysis.tools.MapFeatures;
import org.campagnelab.dl.varanalysis.tools.MapFeaturesArguments;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.dataset.api.MultiDataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
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
    private int index;
    private MultiDataSetPreProcessor preProcessor;

    public MultiDatasetMappedFeaturesIterator(String basename) {
        try {
            Properties cfProperties = new Properties();
            cfProperties.load(new FileReader(basename + ".cfp"));

            miniBatchSize = Integer.parseInt(cfProperties.getProperty("miniBatchSize", "0"));
            numExamples = Integer.parseInt(cfProperties.getProperty("numRecords", "0"));
            inputStream = new FastBufferedInputStream(new FileInputStream(basename + ".cf"));
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
    public void setPreProcessor(MultiDataSetPreProcessor preProcessor)  {
        this.preProcessor = preProcessor;
    }

    @Override
    public boolean resetSupported() {
        return false;
    }


    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        index = 0;
        try {
            inputStream.position(0);
        } catch (IOException e) {
            LOG.error("Unable to reset iterator to position 0");
        }
    }



    @Override
    public boolean hasNext() {
        return index < numExamples;
    }

    byte[] length = new byte[4];
    ByteArrayList content = new ByteArrayList();

    @Override
    public MultiDataSet next() {
        try {
            inputStream.read(length, 0, 4);
        } catch (IOException e) {
            LOG.error("Unable to read length from stream.", e);
        }
        int l = (length[0] << 8 * 3 & 0xFF000000)|
                (length[1] << 8 * 2 & 0x00FF0000) |
                (length[2] << 8     & 0x0000FF00) |
                (length[3]          & 0x000000FF);
        // resize the buffer if needed:
        content.size(l);
        final byte[] elements = content.elements();
        try {
            inputStream.read(elements);
        } catch (IOException e) {
            LOG.error("Unable to read content from stream at index " + index, e);
        }
        MultiDataSet ds = new org.nd4j.linalg.dataset.MultiDataSet();
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
}
