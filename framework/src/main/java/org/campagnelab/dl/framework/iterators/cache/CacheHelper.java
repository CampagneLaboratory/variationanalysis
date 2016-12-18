package org.campagnelab.dl.framework.iterators.cache;

import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.iterators.MultiDatasetMappedFeaturesIterator;
import org.campagnelab.dl.framework.tools.MapMultiDatasetFeatures;
import org.campagnelab.dl.framework.tools.MapMultiDatasetFeaturesArguments;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;

/**
 * A concat iterator that transparently creates a disk cache of the content of the input iterables.
 */
public class CacheHelper<RecordType> {


    /**
     * Return a cached version of the iterator. Either returns a pre-cached iterator, or chaches the iterator
     * and returns the cached version.
     *
     * @param domainDescriptor
     * @param adapter
     * @param cacheName
     * @param cacheN
     * @return A cached iterator.
     */
    public MultiDataSetIterator cache(final DomainDescriptor domainDescriptor,
                                      MultiDataSetIteratorAdapter adapter, String cacheName, int cacheN, int minibatchSize) {


        // determine if cache exists. If it does, use it.

        if (!cacheExists(cacheName, cacheN, true)) {
            // Cache does not exist, we first build it:
            MapMultiDatasetFeatures tool = new MapMultiDatasetFeatures() {
                @Override
                protected DomainDescriptor domainDescriptor() {
                    return domainDescriptor;
                }
            };
            MapMultiDatasetFeaturesArguments arguments = new MapMultiDatasetFeaturesArguments<>();

            arguments.adapter = adapter;
            arguments.outputBasename = cacheName;
            arguments.cacheN = cacheN;
            arguments.domainDescriptor=domainDescriptor;
            arguments.miniBatchSize=minibatchSize;
            tool.setArguments(arguments);
            tool.execute();
        }
        assert cacheExists(cacheName, cacheN, true) : "A cache must exist at this point.";
        return new MultiDatasetMappedFeaturesIterator(cacheName,cacheN);
    }


    /**
     * Check the that cache exists and has the same number of records that indicated in the parameter.
     *
     * @param cacheName
     * @param cacheN
     * @param multiDataSet True when the cache must be MultiDataSet
     * @return
     */
    public static boolean cacheExists(String cacheName, int cacheN, boolean multiDataSet) {
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
            Object descriptor = cfp.getProperty("domainDescriptor");
            if (multiDataSet) {
                if ( descriptor == null) return false;
            }else {
                if ( descriptor != null) return false;
            }

            long cacheNSaved = Long.parseLong(n.toString());
            return (cacheNSaved >= cacheN || cacheN == Integer.MAX_VALUE);
        } catch (FileNotFoundException e) {
            return false;
        } catch (IOException e) {
            return false;
        }
    }


}