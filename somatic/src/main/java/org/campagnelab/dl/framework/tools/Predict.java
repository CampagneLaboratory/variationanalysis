package org.campagnelab.dl.framework.tools;

import com.google.common.collect.Iterables;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.DomainDescriptorLoader;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.iterators.cache.CacheHelper;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

/**
 * A generic Predict tool. Sub-class this abstract class and define a few methods in order to make predictions and
 * evaluation on test sets. See PredictS for an example of sub-class.
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 10/23/16.
 */
public abstract class Predict<RecordType> extends AbstractTool<PredictArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Predict.class);


    public PredictArguments args() {
        return arguments;
    }

    @Override
    public PredictArguments createArguments() {
        return new PredictArguments();
    }

    @Override
    public void execute() {
        PrintWriter resultWriter = null;
        try {
            File modelPath = new File(args().modelPath);
            String modelName = modelPath.getName();

            if (args().toFile) {
                String resultPath = "tests/" + modelName + "/";
                File dir = new File(resultPath);
                // attempt to create the directory here
                dir.mkdirs();

                String resultFilename = resultPath + args().modelName + "-" + args().type + ".tsv";
                System.out.println("Writing predictions to " + resultFilename);
                resultWriter = new PrintWriter(resultFilename, "UTF-8");
            } else {
                resultWriter = new PrintWriter(System.out);
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to create result writer", e);
        }
        try {
            printPredictions(args().modelName, args().modelPath, args().testSet, resultWriter);
        } catch (IOException e) {
            throw new RuntimeException("Unable to perform predictions", e);
        }


    }

    protected DomainDescriptor<RecordType> domainDescriptor;

    private void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                  PrintWriter resutsWriter) throws IOException {


        ModelLoader modelLoader = new ModelLoader(modelPath);
        // we scale features using statistics observed on the training set:
        FeatureMapper featureMapper = modelLoader.loadFeatureMapper(modelLoader.getModelProperties());
        boolean isTrio;
        if (featureMapper.getClass().getCanonicalName().contains("Trio")) {
            //we have a trio mapper, need to output features for a third sample
            isTrio = true;
            System.out.println("setting output to trio mode");
        }
        Model model = modelLoader.loadModel(prefix);
        if (model == null) {
            System.err.println("Cannot load model with prefix: " + prefix);
            System.exit(1);
        }
        domainDescriptor = DomainDescriptorLoader.load(modelPath);
        ProgressLogger pgReadWrite = new ProgressLogger(LOG);
        pgReadWrite.itemsName = prefix;
        final long totalRecords = domainDescriptor.getNumRecords(new String[]{args().testSet});
        pgReadWrite.expectedUpdates = Math.min(args().scoreN,
                totalRecords);
        pgReadWrite.displayFreeMemory = true;
        pgReadWrite.start();

        PredictWithModel<RecordType> predictor = new PredictWithModel<>(domainDescriptor);

        Iterable<RecordType> apply = domainDescriptor.getRecordIterable().apply(evaluationDataFilename);
        Iterable<RecordType> it = Iterables.limit(apply, args().scoreN);
        initializeStats(prefix);
        writeHeader(resutsWriter);

        predictor.makePredictions(it.iterator(),
                model,
                predictionList -> {
                    processPredictions(resutsWriter, predictionList);
                    pgReadWrite.lightUpdate();
                },
                /* stop if */ nProcessed -> nProcessed > args().scoreN
        );

        resutsWriter.close();
        pgReadWrite.stop();
        reportStatistics(prefix);
        modelLoader.writeTestCount(totalRecords);
    }

    /**
     * This method is called after the test set has been observed and statistics evaluated. It prints statistics
     * to the console for the operators to read.
     *
     * @param prefix The model prefix/label (e.g., bestAUC).
     */
    protected abstract void reportStatistics(String prefix);

    /**
     * This method is called for each record of the test set that has been predicted.
     *
     * @param predictionList The list of predictions made for a test record. Typically one prediction for each model output.
     * @param resutsWriter   Where predictions can be written in tab delited format.
     */
    protected abstract void processPredictions(PrintWriter resutsWriter, List<Prediction> predictionList);

    /**
     * This method is called when we need to write the header to the results.
     *
     * @param resutsWriter
     */
    protected abstract void writeHeader(PrintWriter resutsWriter);

    /**
     * This method must allocate any statistic caculator necessary to evaluate performance on the test set and store
     * these instances in fields of the sub-class.
     *
     * @param prefix The model prefix/label (e.g., bestAUC).
     */
    protected abstract void initializeStats(String prefix);


}