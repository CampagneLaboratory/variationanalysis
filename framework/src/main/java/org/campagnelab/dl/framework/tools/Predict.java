package org.campagnelab.dl.framework.tools;

import com.google.common.collect.Iterables;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.domains.DomainDescriptorLoader;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.ConditionRecordingTool;
import org.deeplearning4j.nn.api.Model;
import org.nd4j.linalg.dataset.api.MultiDataSet;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A generic Predict tool. Sub-class this abstract class and define a few methods in order to make predictions and
 * evaluation on test sets. See PredictS for an example of sub-class.
 *
 * @author Fabien Campagne
 *         Created by fac2003 on 10/23/16.
 */
public abstract class Predict<RecordType> extends ConditionRecordingTool<PredictArguments> {

    static private Logger LOG = LoggerFactory.getLogger(Predict.class);
    protected String modelTime;
    protected String modelPrefix;
    protected String testSetBasename;

    public PredictArguments args() {
        return arguments;
    }

    @Override
    public PredictArguments createArguments() {
        return new PredictArguments();
    }

    @Override
    public void execute() {
        if (args().deviceIndex != null) {
            Nd4j.getAffinityManager().attachThreadToDevice(Thread.currentThread(), args().deviceIndex);
        }
        PrintWriter resultWriter;
        PrintWriter outputWriter;
        boolean outputFileExists;
        try {
            File modelPath = new File(args().modelPath);
            String modelTime = modelPath.getName();
            File file = new File(args().outputFile);
            outputFileExists = file.exists() && file.length() > 0;
            if (args().toFile) {
                String resultPath = "predictions";
                File dir = new File(resultPath);
                // attempt to create the directory here
                dir.mkdirs();
                String testSetName = FilenameUtils.getBaseName(args().testSet);
                String resultFilename = String.format("%s/%s-%s-%s-%s.tsv", resultPath, modelTime, args().modelName,
                        args().type, testSetName);
                System.out.println("Writing predictions to " + resultFilename);
                resultWriter = new PrintWriter(resultFilename, "UTF-8");
                outputWriter = new PrintWriter(new FileWriter(args().outputFile, true));
                this.modelTime = modelTime;
                this.modelPrefix = args().modelName;
                this.testSetBasename = testSetName;
            } else {
                resultWriter = new PrintWriter(System.out);
                outputWriter = new PrintWriter(System.out);
                outputFileExists = false;
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to create result writer", e);
        }
        try {
            printPredictions(args().modelName, args().modelPath, args().testSet, resultWriter, outputWriter,
                    outputFileExists);
        } catch (IOException e) {
            throw new RuntimeException("Unable to perform predictions", e);
        }


    }

    protected DomainDescriptor<RecordType> domainDescriptor;

    private void printPredictions(String prefix, String modelPath, String evaluationDataFilename,
                                  PrintWriter resutsWriter, PrintWriter outputWriter,
                                  boolean outputFileExists) throws IOException {


        ModelLoader modelLoader = new ModelLoader(modelPath);
        String modelTag = modelLoader.getModelProperties().getProperty("tag");
        if (!outputFileExists) {
            outputWriter.append("tag\tprefix");
            for (String metricName : createOutputHeader()) {
                outputWriter.append(String.format("\t%s", metricName));
            }
            outputWriter.append("\targuments\n");
        }
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
        pgReadWrite.itemsName = "sites";
        final long totalRecords = domainDescriptor.getNumRecords(new String[]{args().testSet});
        pgReadWrite.expectedUpdates = Math.min(args().scoreN,
                totalRecords);
        pgReadWrite.displayFreeMemory = true;


        PredictWithModel<RecordType> predictor = new PredictWithModel<RecordType>(domainDescriptor);

        Iterable<RecordType> apply = domainDescriptor.getRecordIterable().apply(evaluationDataFilename);
        Iterable<RecordType> itAdapter = Iterables.limit(apply, args().scoreN);
        Iterable<RecordType> recordsIterable = Iterables.limit(domainDescriptor.getRecordIterable().apply(evaluationDataFilename), args().scoreN);
        initializeStats(prefix);
        writeHeader(resutsWriter);
        final int miniBatchSize = args().miniBatchSize;
        MultiDataSetIteratorAdapter<RecordType> adapter = new MultiDataSetIteratorAdapter<RecordType>(itAdapter,
                miniBatchSize, domainDescriptor, false, null) {
            @Override
            public String getBasename() {
                return FilenameUtils.getBaseName(args().testSet);
            }
        };
        List<RecordType> records = new ArrayList<RecordType>();
        Iterator<RecordType> recordIterator = recordsIterable.iterator();
        int index = 0;
        int adapterIndex=0;
        pgReadWrite.start();
        while (adapter.hasNext() && recordIterator.hasNext()) {

            MultiDataSet dataset = adapter.next();
            final int datasetSize = dataset.getFeatures(0).rows();
            adapterIndex++;
            records.clear();
            for (int exampleIndex = 0; exampleIndex < datasetSize; exampleIndex++) {
                if (!recordIterator.hasNext()) {
                    break;
                }
                records.add(recordIterator.next());
            }

            if (records.size()==datasetSize) {
                index = predictor.makePredictions(dataset,
                        records, model,
                        recordPredictions -> {
                            processPredictions(resutsWriter, recordPredictions.record,
                                    recordPredictions.predictions);
                            pgReadWrite.update(records.size());
                        },
                /* stop if */ nProcessed -> nProcessed > args().scoreN, index
                );
            } else{
                System.out.printf("dataset #examples %d and # records (%d) must match. Unable to obtain records for some examples in minibatch. Aborting. ",
                        datasetSize, records.size());
                break;
            }

        }


        resutsWriter.close();
        outputWriter.append(String.format("%s\t%s", modelTag, prefix));
        for (double metric : createOutputStatistics()) {
            outputWriter.append(String.format("\t%f", metric));
        }
        outputWriter.append("\t" + getAllCommandLineArguments());
        outputWriter.append("\n");
        outputWriter.close();
        pgReadWrite.stop();
        reportStatistics(prefix);
        System.out.println("Model: " + modelPath + " tag:" + modelTag);
        modelLoader.writeTestCount(totalRecords);
    }

    /**
     * This method is called after the test set has been observed and statistics evaluated via processPredictions.
     * It sets statistics on the whole test set, which are then written tab-delimited to a file.
     *
     * @return array of statistics to set; should be in same order as outputHeader
     */
    protected abstract double[] createOutputStatistics();

    /**
     * This method is called only if an output file for statistics has not been created yet. It sets the fields in
     * the header file, to then be written to the output file.
     *
     * @return array of metric names to set
     */
    protected abstract String[] createOutputHeader();

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
    protected abstract void processPredictions(PrintWriter resutsWriter, RecordType record, List<Prediction> predictionList);

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