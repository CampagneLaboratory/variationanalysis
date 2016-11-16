package org.campagnelab.dl.varanalysis.learning.performance;

import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.varanalysis.mappers.FeatureMapper;
import org.campagnelab.dl.varanalysis.mappers.LabelMapper;
import org.campagnelab.dl.varanalysis.learning.iterators.BaseInformationIterator;
import org.campagnelab.dl.varanalysis.learning.iterators.FirstNIterator;
import org.campagnelab.dl.varanalysis.stats.AUCHelper;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.dataset.api.iterator.CachingDataSetIterator;
import org.nd4j.linalg.dataset.api.iterator.cache.InMemoryDataSetCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * Helper class to estimate performance on a test set, or part of it.
 * Created by fac2003 on 7/15/16.
 */
public class MeasurePerformance {
    private int aucClipMaxObservations;
    private int scoreN = Integer.MAX_VALUE;
    private AUCHelper helper;
    static private Logger LOG = LoggerFactory.getLogger(MeasurePerformance.class);

    public MeasurePerformance(int scoreN, String datafilePath, int miniBatchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {

        this.aucClipMaxObservations = Integer.MAX_VALUE;
        this.scoreN = scoreN;
        ProgressLogger pg = new ProgressLogger(LOG);
        cachedIterator = new CachingDataSetIterator(
                new FirstNIterator(new BaseInformationIterator(datafilePath, miniBatchSize, featureMapper, labelMapper), scoreN),
                new InMemoryDataSetCache());
        pg.expectedUpdates = cachedIterator.numExamples()/miniBatchSize;
        pg.itemsName = "validation records";
        pg.start();
        while (cachedIterator.hasNext()) {
            cachedIterator.next();
            pg.lightUpdate();
        }
        pg.stop();
        cachedIterator.reset();
        helper = new AUCHelper();
    }

    public MeasurePerformance(int scoreN, int aucClipMaxObservations, String datafilePath, int miniBatchSize, FeatureMapper featureMapper, LabelMapper labelMapper) throws IOException {
        this(scoreN, datafilePath, miniBatchSize, featureMapper, labelMapper);
        this.aucClipMaxObservations = aucClipMaxObservations;
    }

    CachingDataSetIterator cachedIterator;


    public double estimateAUC(MultiLayerNetwork model) throws IOException {

        cachedIterator.reset();
        return helper.estimate(cachedIterator, model, this.aucClipMaxObservations,
                prediction -> {
                },
                numScored -> numScored > scoreN);
    }


}
