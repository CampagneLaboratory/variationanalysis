package org.campagnelab.dl.varanalysis.learning.iterators;

import org.campagnelab.dl.model.utils.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.tools.MapFeatures;
import org.campagnelab.dl.varanalysis.tools.MapFeaturesArguments;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

/**
 * A concat iterator that transparently creates a disk cache of the content of the input iterators.
 */
public class CachedConcatIterator<T> implements NamedDataSetIterator {


    MappedFeaturesIterator delegate;
    private List<NamedDataSetIterator> iterators;

    public CachedConcatIterator(NamedDataSetIterator iterator, int minbatchSize, FeatureMapper featureMapper, long cacheN) {
        this(Collections.singletonList(iterator), minbatchSize, featureMapper, cacheN);
    }

    public CachedConcatIterator(List<NamedDataSetIterator> iterators, int minbatchSize, FeatureMapper<T> featureMapper, long cacheN) {
        this.iterators = iterators;
        // determine if cache exists. If it does, use it.
        String cacheName = buildCacheName(iterators);
        if (!cacheExists(cacheName, cacheN)) {
            // Cache does not exist, we first build it:
            MapFeatures tool = new MapFeatures();
            MapFeaturesArguments arguments = new MapFeaturesArguments();

            for (NamedDataSetIterator it : iterators) {
                arguments.trainingSets.add(it.getBasename());
            }
            arguments.miniBatchSize = minbatchSize;
            arguments.featureMapperClassname = featureMapper.getClass().getCanonicalName();
            arguments.outputBasename = cacheName;
            arguments.cacheN = cacheN;
            tool.setArguments(arguments);
            tool.execute();
        }
        assert cacheExists(cacheName, cacheN) : "A cache must exist at this point.";
        delegate = new MappedFeaturesIterator(cacheName);
    }


    @Override
    public DataSet next(int miniBatchSize) {
        return delegate.next(miniBatchSize);
    }

    @Override
    public int totalExamples() {
        return delegate.totalExamples();
    }

    @Override
    public int inputColumns() {
        return delegate.inputColumns();
    }

    @Override
    public int totalOutcomes() {
        return delegate.totalOutcomes();
    }

    @Override
    public boolean resetSupported() {
        return delegate.resetSupported();
    }

    public boolean asyncSupported() {
        return false;
    }

    @Override
    public void reset() {
        delegate.reset();
    }

    @Override
    public int batch() {
        return delegate.batch();
    }

    @Override
    public int cursor() {
        return delegate.cursor();
    }

    @Override
    public int numExamples() {
        return delegate.numExamples();
    }

    @Override
    public void setPreProcessor(DataSetPreProcessor dataSetPreProcessor) {
        delegate.setPreProcessor(dataSetPreProcessor);
    }

    @Override
    public DataSetPreProcessor getPreProcessor() {
        return delegate.getPreProcessor();
    }

    @Override
    public List<String> getLabels() {
        return delegate.getLabels();
    }

    @Override
    public boolean hasNext() {
        return delegate.hasNext();
    }

    @Override
    public DataSet next() {
        return delegate.next();
    }

    /**
     * Check the that cache exists and has the same number of records that indicated in the parameter.
     *
     * @param cacheName
     * @param cacheN
     * @return
     */
    private boolean cacheExists(String cacheName, long cacheN) {
        boolean cacheExists = new File(cacheName + ".cf").exists() & new File(cacheName + ".cfp").exists();
        if (!cacheExists) {
            return false;
        }
        try {
            // check that the number of cached items matches the value of cacheN on the command line:
            Properties cfp = new Properties();
            cfp.load(new FileReader(new File(cacheName + ".cfp")));
            Object n = cfp.getProperty("numRecords");
            if (n == null) return false;
            long cacheNSaved = Long.parseLong(n.toString());
            return (cacheNSaved >= cacheN || cacheN == Long.MAX_VALUE);
        } catch (FileNotFoundException e) {
            return false;
        } catch (IOException e) {
            return false;
        }
    }

    private String buildCacheName(List<NamedDataSetIterator> iterators) {
        String cacheName;// only one input, use its name as cache name:
        if (iterators.size() == 1) {
            cacheName = iterators.get(0).getBasename();
        } else {
            long hashcode = 8723872838723L;
            for (NamedDataSetIterator it : iterators) {
                hashcode ^= it.getBasename().hashCode();
            }
            cacheName = "multiset-" + Long.toString(hashcode);

        }
        return cacheName;
    }

    @Override
    public String getBasename() {
        return buildCacheName(iterators);
    }
}