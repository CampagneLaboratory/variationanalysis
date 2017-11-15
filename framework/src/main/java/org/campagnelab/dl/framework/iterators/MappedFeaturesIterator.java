package org.campagnelab.dl.framework.iterators;

import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.io.FastBufferedInputStream;
import org.apache.commons.lang.NotImplementedException;
import org.campagnelab.dl.framework.mappers.LabelMapper;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;
import org.nd4j.linalg.dataset.api.iterator.DataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Properties;

/**
 * An Iterator over mapped features (.cf/.cfp files).
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 11/2/16.
 */
public class MappedFeaturesIterator implements DataSetIterator {
    static private Logger LOG = LoggerFactory.getLogger(MappedFeaturesIterator.class);

    private final int numExamples;
    private final int miniBatchSize;
    private final int numColumns;
    private final LabelMapper labelMapper;
    private final FastBufferedInputStream inputStream;
    private int index;
    private DataSetPreProcessor preProcessor;

    public MappedFeaturesIterator(String basename) {
        try {
            Properties cfProperties = new Properties();
            cfProperties.load(new FileReader(basename + ".cfp"));
            miniBatchSize = Integer.parseInt(cfProperties.getProperty("miniBatchSize", "0"));
            numExamples = Integer.parseInt(cfProperties.getProperty("numRecords", "0"));
            numColumns = Integer.parseInt(cfProperties.getProperty("numFeatures", "0"));
            String labelMapperClassname = cfProperties.getProperty("labelMapper");
            labelMapper = (LabelMapper) Class.forName(labelMapperClassname).newInstance();
            inputStream = new FastBufferedInputStream(new FileInputStream(basename + ".cf"));
        } catch (Exception e) {
            throw new RuntimeException("Unable to create MappedFeaturesIterator ", e);
        }
    }

    @Override
    public DataSet next(int miniBatchSize) {
        if (miniBatchSize != this.miniBatchSize) {
            throw new IllegalArgumentException("numExamples must match the cached minibatchSize: " + miniBatchSize);
        }
        return next();
    }

    @Override
    public int totalExamples() {
        return numExamples;
    }

    @Override
    public int inputColumns() {
        return numColumns;
    }

    @Override
    public int totalOutcomes() {
        return labelMapper.numberOfLabels();
    }

    @Override
    public boolean resetSupported() {
        return true;
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
    public int batch() {
        return miniBatchSize;
    }

    @Override
    public int cursor() {
        return index;
    }

    @Override
    public int numExamples() {
        return numExamples;
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor dataSetPreProcessor) {
        this.preProcessor = dataSetPreProcessor;
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
        return preProcessor;
    }

    @Override
    public List<String> getLabels() {
        throw new NotImplementedException();
    }

    @Override
    public boolean hasNext() {
        return index < numExamples;
    }

    byte[] length = new byte[4];
    ByteArrayList content = new ByteArrayList();

    @Override
    public DataSet next() {
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
        DataSet ds = new DataSet();
        try (ByteArrayInputStream from = new ByteArrayInputStream(elements)) {
            ds.load(from);
        } catch (IOException e) {
            LOG.error("Unable to load dataset at index " + index, e);
        }
        if (preProcessor != null) {
            preProcessor.preProcess(ds);
        }
        index += ds.numExamples();
        return ds;

    }
}
