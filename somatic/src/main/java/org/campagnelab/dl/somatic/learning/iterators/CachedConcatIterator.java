package org.campagnelab.dl.somatic.learning.iterators;

import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.somatic.tools.MapFeatures;
import org.campagnelab.dl.somatic.tools.MapFeaturesArguments;
import org.nd4j.linalg.dataset.DataSet;
import org.nd4j.linalg.dataset.api.DataSetPreProcessor;

import java.util.Collections;
import java.util.List;

/**
 * A concat iterator that transparently creates a disk cache of the content of the input iterables.
 */
public class CachedConcatIterator<T> implements NamedDataSetIterator {


    MappedFeaturesIterator delegate;
    private List<NamedDataSetIterator> iterators;

    public CachedConcatIterator(NamedDataSetIterator iterator, int minbatchSize,
                                FeatureMapper featureMapper, int cacheN) {
        this(Collections.singletonList(iterator), minbatchSize, featureMapper, cacheN);
    }

    public CachedConcatIterator(List<NamedDataSetIterator> iterators, int minbatchSize,
                                FeatureMapper<T> featureMapper, int cacheN) {
        this.iterators = iterators;
        // determine if cache exists. If it does, use it.
        String cacheName = buildCacheName(iterators);
        if (!CacheHelper.cacheExists(cacheName, cacheN,false)) {
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
        assert CacheHelper.cacheExists(cacheName, cacheN,false) : "A cache must exist at this point.";
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