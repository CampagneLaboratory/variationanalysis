package org.campagnelab.dl.framework.tools;

import com.google.common.collect.Iterables;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.floats.FloatArraySet;
import it.unimi.dsi.fastutil.floats.FloatSet;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.architecture.graphs.ComputationGraphAssembler;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.gpu.InitializeGpu;
import org.campagnelab.dl.framework.gpu.ParameterPrecision;
import org.campagnelab.dl.framework.iterators.MultiDataSetIteratorAdapter;
import org.campagnelab.dl.framework.iterators.cache.CacheHelper;
import org.campagnelab.dl.framework.iterators.cache.FullyInMemoryCache;
import org.campagnelab.dl.framework.mappers.FeatureMapper;
import org.campagnelab.dl.framework.models.ComputationGraphSaver;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.campagnelab.dl.framework.performance.Metric;
import org.campagnelab.dl.framework.performance.PerformanceLogger;
import org.campagnelab.dl.framework.performance.PerformanceMetricDescriptor;
import org.campagnelab.dl.framework.tools.arguments.ConditionRecordingTool;
import org.campagnelab.dl.framework.training.ParallelTrainerOnGPU;
import org.campagnelab.dl.framework.training.SequentialTrainer;
import org.campagnelab.dl.framework.training.Trainer;
import org.deeplearning4j.api.storage.StatsStorage;
import org.deeplearning4j.earlystopping.EarlyStoppingResult;
import org.deeplearning4j.nn.api.Layer;
import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.deeplearning4j.ui.api.UIServer;
import org.deeplearning4j.ui.stats.StatsListener;
import org.deeplearning4j.ui.storage.FileStatsStorage;
import org.deeplearning4j.ui.storage.InMemoryStatsStorage;
import org.nd4j.linalg.api.buffer.DataBuffer;
import org.nd4j.linalg.api.buffer.util.DataTypeUtil;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.dataset.api.iterator.MultiDataSetIterator;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * An abstract tool to train computational graphs. Implements early stopping. This class defines
 * several abstract methods that must be implemented to adapt training to different problems.
 */
public abstract class TrainModel<RecordType> extends ConditionRecordingTool<TrainingArguments> {

    static private Logger LOG = LoggerFactory.getLogger(TrainModel.class);

    private String directory;
    private double bestScore;
    private long time;

    protected DomainDescriptor<RecordType> domainDescriptor;
    private String bestMetricName;

    protected abstract DomainDescriptor<RecordType> domainDescriptor();

    protected PerformanceLogger performanceLogger;

    protected FeatureMapper featureMapper = null;
    private ComputationGraph computationGraph;
    private CacheHelper<RecordType> cacheHelper = new CacheHelper<>();


    @Override
    public void execute() {
        InitializeGpu.initialize();
        if (args().getTrainingSets().length == 0) {
            System.err.println("You must provide training datasets.");
        }
        domainDescriptor = domainDescriptor();
        try {
            featureMapper = domainDescriptor().getFeatureMapper("input");

            execute(featureMapper, args().getTrainingSets(), args().miniBatchSize);
        } catch (IOException e) {
            System.err.println("An exception occured. Details may be provided below");
            e.printStackTrace();
        }

    }

    public void execute(FeatureMapper featureCalculator, String trainingDataset[], int miniBatchSize) throws IOException {

        if (args().deviceIndex != null && !args().parallel) {
            Nd4j.getAffinityManager().attachThreadToDevice(Thread.currentThread(), args().deviceIndex);
        }
        if (args().previousModelPath != null) {
            System.out.println(String.format("Resuming training with %s model parameters from %s %n", args().previousModelName, args().previousModelPath));
        }
        if ("FP16".equals(args().precision)) {
            precision = ParameterPrecision.FP16;
            System.out.println("Parameter precision set to FP16.");
        }
        if ("FP16".equals(args().precision)) {
            DataTypeUtil.setDTypeForContext(DataBuffer.Type.HALF);
        }
        time = new Date().getTime();

        System.out.println("epochs: " + args().maxEpochs);
        System.out.println(featureCalculator.getClass().getTypeName());
        directory = "models/" + Long.toString(time);
        FileUtils.forceMkdir(new File(directory));
        System.out.println("model directory: " + new File(directory).getAbsolutePath());

        // Assemble the computational graph:
        bestMetricName = "best" + domainDescriptor.performanceDescritor().earlyStoppingMetric();

        ComputationGraphAssembler assembler = domainDescriptor.getComputationalGraph();
        assert assembler != null : "Computational Graph assembler must be defined.";
        assembler.setArguments(args());
        Map<String, Boolean> padded = new HashMap<>();
        for (String inputName : assembler.getInputNames()) {
            int[] domainDescriptorNumInputs = domainDescriptor.getNumInputs(inputName).clone();
            boolean domainPadEos = (args().previousModelPretraining) &&
                    ((args().eosIndex != null && args().eosIndex == domainDescriptorNumInputs[0])
                            || args().eosIndex == null);
            padded.put(inputName, domainPadEos);
            if (domainPadEos) domainDescriptorNumInputs[0]++;
            assembler.setNumInputs(inputName, domainDescriptorNumInputs);
        }
        domainDescriptor.setInputsPaddedEos(padded);
        for (String outputName : assembler.getOutputNames()) {
            assembler.setNumOutputs(outputName, domainDescriptor.getNumOutputs(outputName));
            assembler.setLossFunction(outputName, domainDescriptor.getOutputLoss(outputName));
        }
        for (String componentName : assembler.getComponentNames()) {
            assembler.setNumHiddenNodes(componentName, domainDescriptor.getNumHiddenNodes(componentName));
        }

        computationGraph = assembler.createComputationalGraph(domainDescriptor);
        computationGraph.init();
        if (args().addUiListener) {
            UIServer uiServer = UIServer.getInstance();
            StatsStorage statsStorage = args().uiStatsFile == null
                    ? new InMemoryStatsStorage()
                    : new FileStatsStorage(new File(args().uiStatsFile));
            uiServer.attach(statsStorage);
            computationGraph.setListeners(new StatsListener(statsStorage));
        }
        if (args().previousModelPath != null) {
            // Load the parameters of a previously trained model and set them on the new model to continue
            // training where we left it off. Note that models must have the same architecture or setting
            // parameters will fail.

            ModelLoader loader = new ModelLoader(args().previousModelPath);
            Model savedNetwork = loader.loadModel(args().previousModelName);
            ComputationGraph savedGraph = savedNetwork instanceof ComputationGraph ?
                    (ComputationGraph) savedNetwork :
                    null;
            if (savedNetwork == null || savedGraph == null || savedGraph.getUpdater() == null || savedGraph.params() == null) {
                System.err.println("Unable to load model or updater from " + args().previousModelPath);
            } else {
                computationGraph.setUpdater(savedGraph.getUpdater());
                computationGraph.setParams(savedNetwork.params());
            }
        }
        //Print the  number of parameters in the graph (and for each layer)
        Layer[] layers = computationGraph.getLayers();
        int totalNumParams = 0;
        for (int i = 0; i < layers.length; i++) {
            int nParams = layers[i].numParams();
            System.out.println("Number of parameters in layer " + i + ": " + nParams);
            totalNumParams += nParams;
        }
        System.out.println("Total number of network parameters: " + totalNumParams);

        writeProperties();
        performanceLogger = new PerformanceLogger(directory);
        PerformanceMetricDescriptor<RecordType> perfDescriptor = domainDescriptor.performanceDescritor();

        Metric[] metrics = new Metric[perfDescriptor.performanceMetrics().length];
        for (int i = 0; i < metrics.length; i++) {
            String metricName = perfDescriptor.performanceMetrics()[i];
            metrics[i] = new Metric(metricName, perfDescriptor.largerValueIsBetterPerformance(metricName));
        }
        performanceLogger.definePerformances(metrics);
        EarlyStoppingResult<ComputationGraph> result = train();

        //Print out the results:
        System.out.println("Termination reason: " + result.getTerminationReason());
        System.out.println("Termination details: " + result.getTerminationDetails());
        System.out.println("Total epochs: " + result.getTotalEpochs());
        System.out.println("Best epoch number: " + result.getBestModelEpoch());

        for (String metricName : perfDescriptor.performanceMetrics()) {
            System.out.println(metricName + " at best epoch: " + performanceLogger.getBest(metricName));
        }
        writeProperties();
        writeBestScoreFile();
        System.out.println("Model completed, saved at time: " + time);
        performanceLogger.write();
        for (String metric : domainDescriptor.performanceDescritor().performanceMetrics()) {
            resultValues().put(metric, performanceLogger.getBest(metric));

        }
        resultValues().put("bestModelEpoch", performanceLogger.getBestEpoch(bestMetricName));
        resultValues().put("model-time", time);
    }


    protected void writeBestScoreFile() throws IOException {

        FileWriter scoreWriter = new FileWriter(directory + "/bestScore");
        scoreWriter.append(Double.toString(performanceLogger.getBestScore()));
        scoreWriter.close();
    }

    protected void writeProperties() throws IOException {
        ModelPropertiesHelper mpHelper = new ModelPropertiesHelper();
        ComputationGraphAssembler assembler = domainDescriptor.getComputationalGraph();
        appendProperties(assembler, mpHelper);
        mpHelper.addProperties(getReaderProperties(args().trainingSets.get(0)));
        mpHelper.put("domainDescriptor", domainDescriptor.getClass().getCanonicalName());
        mpHelper.put("tag", getTag());
        domainDescriptor().putProperties(mpHelper.getProperties());
        mpHelper.writeProperties(directory);
        domainDescriptor.writeProperties(directory);
    }


    protected static int numLabels(INDArray labels) {
        FloatSet set = new FloatArraySet();
        for (int i = 0; i < labels.size(0); i++) {
            set.add(labels.getFloat(i));
        }
        return set.size();
    }


    public void appendProperties(ComputationGraphAssembler assembler, ModelPropertiesHelper helper) {
        // give a chance to the assembler to save information to the model properties to describe the architecture
        // used for training:
        assembler.saveProperties(helper);

        //save the rest of the arguments:
        helper.setFeatureCalculator(featureMapper);
        helper.setLearningRate(args().learningRate);
        helper.setDropoutRate(args().dropoutRate);
        helper.setMiniBatchSize(args().miniBatchSize);
        // mpHelper.setBestScore(bestScore);
        helper.setNumEpochs(args().maxEpochs);
        helper.setNumTrainingSets(args().trainingSets.size());
        helper.setTime(time);
        helper.setSeed(args().seed);

        helper.setEarlyStopCriterion(args().stopWhenEpochsWithoutImprovement);
        helper.setRegularization(args().regularizationRate);
        helper.setPrecision(precision);
        helper.put("allArguments", getAllCommandLineArguments());
    }

    ParameterPrecision precision = ParameterPrecision.FP32;

    public abstract Properties getReaderProperties(String trainingSet) throws IOException;


    protected EarlyStoppingResult<ComputationGraph> train() throws IOException {
        String validationDatasetFilename = args().validationSet;
        //check validation file for error
        if (!(new File(validationDatasetFilename).exists())) {
            throw new IOException("Validation file not found! " + validationDatasetFilename);
        }
        //Do training, and then generate and print samples from network
        int miniBatchNumber = 0;
        boolean init = true;
        bestScore = Double.MAX_VALUE;
        ComputationGraphSaver saver = new ComputationGraphSaver(directory);
        int iter = 0;
        Map<Integer, Double> scoreMap = new HashMap<Integer, Double>();
        System.out.println("errorEnrichment=" + args().errorEnrichment);

        performanceLogger.setCondition(args().experimentalCondition);
        long numExamplesUsed = 0;
        int notImproved = 0;

        System.out.flush();
        PerformanceMetricDescriptor perfDescriptor = domainDescriptor.performanceDescritor();
        String validationMetricName = perfDescriptor.earlyStoppingMetric();
        double bestValue = initializePerformance(perfDescriptor, validationMetricName);
        int epoch;

        // Assemble the training iterator from the concatenation of individual training set iterables:
        Iterable<RecordType> inputIterable = Iterables.concat(
                args().trainingSets.stream().map(
                        filename -> domainDescriptor.getRecordIterable().apply(filename)).collect(
                        Collectors.toList()));
        Iterable<RecordType> recordIterable = Iterables.limit(inputIterable, args().numTraining);
        final int miniBatchSize = args().miniBatchSize;
        MultiDataSetIteratorAdapter<RecordType> adapter = new MultiDataSetIteratorAdapter<RecordType>(recordIterable,
                miniBatchSize, domainDescriptor, args().previousModelPretraining, args().eosIndex) {
            @Override
            public String getBasename() {
                return buildBaseName(args().trainingSets);
            }
        };

        boolean useCache = !args().ignoreCache;
        MultiDataSetIterator iterator = useCache ? cacheHelper.cache(domainDescriptor,
                adapter, adapter.getBasename(),
                args().numTraining, args().miniBatchSize) :
                adapter;
        if (args().memoryCacheTraining()) {
            iterator = new FullyInMemoryCache(iterator);
            // force loading immediately:
            LOG.warn("Loading training set in memory.");
            iterator.reset();
            LOG.warn("Done.");
        }
        // MultiDataSetIterator iterator=adapter;
        final long numRecords = Math.min(args().numTraining, domainDescriptor.getNumRecords(args().getTrainingSets()));
        int miniBatchesPerEpoch = (int) (numRecords / args().miniBatchSize);
        System.out.printf("Training with %d minibatches per epoch%n", miniBatchesPerEpoch);
        MultiDataSetIterator validationIterator = readValidationSet();
        System.out.println("Finished loading validation records.");


        if (args().buildCacheAndStop) {
            System.out.println("Cache has been built. Exiting now since --build-cache-then-stop was used.");
            System.exit(0);
        }
        ProgressLogger pgEpoch = new ProgressLogger(LOG);
        pgEpoch.displayLocalSpeed = true;
        pgEpoch.itemsName = "epoch";
        pgEpoch.expectedUpdates = args().maxEpochs;

        switch (args().trackingStyle) {
            case PERFS:
                System.out.println(performanceLogger.getMetricHeader());
                break;
            case SPEED:
                pgEpoch.start();
                break;
            default:
                System.out.println("Unsupported tracking style: " + args().trackingStyle);
        }

        Trainer trainer = args().parallel ? new ParallelTrainerOnGPU(computationGraph, args().miniBatchSize,
                (int) domainDescriptor.getNumRecords(args().getTrainingSets())) :
                new SequentialTrainer();
        trainer.setLogSpeed(args().trackingStyle == TrainingArguments.TrackStyle.SPEED);
        for (epoch = 0; epoch < args().maxEpochs; epoch++) {
            ProgressLogger pg = new ProgressLogger(LOG);
            pg.itemsName = "mini-batch";
            iter = 0;
            pg.expectedUpdates = miniBatchesPerEpoch; // one iteration processes miniBatchIterator elements.
            if (args().trackingStyle == TrainingArguments.TrackStyle.SPEED) {
                pg.start();
            }
            // train the graph with the content of the iterator:
            numExamplesUsed += trainer.train(computationGraph, iterator, pg);

            //save latest after the end of an epoch:
            double trainingScore = computationGraph.score();
            saver.saveLatestModel(computationGraph, trainingScore);
            writeProperties();
            writeBestScoreFile();
            if (epoch % args().validateEvery == 0) {

                // estimate all performance metrics. Note that we do a pass over the validation set for each metric:
                // (avoid using several metrics).
                double validationMetricValue = initializePerformance(perfDescriptor, perfDescriptor.earlyStoppingMetric());
                DoubleArrayList metricValues = new DoubleArrayList();


                validationIterator.reset();
                assert validationIterator.hasNext() : "validation iterator must have datasets. Make sure the latest release of Goby is installed in the maven repo.";
                final double[] performanceValues = perfDescriptor.estimateMetric(computationGraph,
                        validationIterator, args().numValidation, perfDescriptor.performanceMetrics());
                metricValues = DoubleArrayList.wrap(performanceValues);

                validationMetricValue = findMetricValue(perfDescriptor.earlyStoppingMetric(),
                        perfDescriptor.performanceMetrics(),
                        performanceValues);

                performanceLogger.logMetrics("epochs", numExamplesUsed, epoch, metricValues.toDoubleArray());
                performanceLogger.logTrainingScore("epochs", epoch, trainingScore);
                if (args().trackingStyle == TrainingArguments.TrackStyle.PERFS) {
                    performanceLogger.show("epochs");
                }
                //System.out.println(metricValues);
                if (!Double.isNaN(bestValue) &&
                        (perfDescriptor.largerValueIsBetterPerformance(validationMetricName) && validationMetricValue > bestValue) ||
                        (!perfDescriptor.largerValueIsBetterPerformance(validationMetricName) && validationMetricValue < bestValue)) {
                    saver.saveModel(computationGraph, "best" + validationMetricName);
                    bestValue = validationMetricValue;

                    performanceLogger.logMetrics(bestMetricName, numExamplesUsed, epoch, metricValues.toDoubleArray());
                    notImproved = 0;
                } else {
                    notImproved++;
                }
                if (notImproved > args().stopWhenEpochsWithoutImprovement) {
                    // we have not improved after earlyStopCondition epoch, time to stop.
                    break;
                }
            }
            if (args().trackingStyle == TrainingArguments.TrackStyle.SPEED) {
                pg.stop();
                pgEpoch.updateAndDisplay();
            }
            iterator.reset();    //Reset iterator for another epoch
            performanceLogger.write();
            //addCustomOption("--error-enrichment", args().errorEnrichment);
            //addCustomOption("--num-errors-added", args().numErrorsAdded);
        }
        pgEpoch.stop();
        return new EarlyStoppingResult<ComputationGraph>(EarlyStoppingResult.TerminationReason.EpochTerminationCondition,
                "not early stopping", scoreMap, performanceLogger.getBestEpoch(bestMetricName), bestScore, args().maxEpochs, computationGraph);
    }


    private double findMetricValue(String lookupName, String[] metricNames, double[] performanceValues) {
        int i = 0;
        for (String name : metricNames) {
            if (lookupName.equals(name)) {
                return performanceValues[i];
            }
            i++;
        }
        throw new RuntimeException("Metric name not found: " + lookupName);
    }

    private String buildBaseName(List<String> trainingSets) {
        String cacheName;// only one input, use its name as cache name:
        if (trainingSets.size() == 1) {

            cacheName = FilenameUtils.removeExtension(trainingSets.get(0));
            ;
        } else {
            long hashcode = 8723872838723L;
            for (String name : trainingSets) {
                hashcode ^= FilenameUtils.getBaseName(name).hashCode();
            }
            cacheName = "multiset-" + Long.toString(hashcode);

        }
        return cacheName;
    }


    private double initializePerformance(PerformanceMetricDescriptor perfDescriptor, String validationMetricName) {
        return perfDescriptor.largerValueIsBetterPerformance(validationMetricName) ?
                Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
    }

    private MultiDataSetIterator readValidationSet() {
        Iterable<RecordType> validationRecords = domainDescriptor.getRecordIterable().apply(args().validationSet);
        try {
            MultiDataSetIteratorAdapter<RecordType> adapter = new MultiDataSetIteratorAdapter<RecordType>(validationRecords,
                    args().miniBatchSize, domainDescriptor, args().previousModelPretraining, args().eosIndex) {
                @Override
                public String getBasename() {
                    return args().validationSet;
                }
            };
            MultiDataSetIterator iterator = args().ignoreCache ? adapter : cacheHelper.cache(domainDescriptor,
                    adapter, adapter.getBasename(),
                    args().numValidation, args().miniBatchSize);
            if (args().memoryCacheValidation()) {
                iterator = new FullyInMemoryCache(iterator);
            }
            return iterator;
        } catch (IOException e) {
            throw new RuntimeException("Unable to load validation records from " + args().validationSet);
        }

    }


}
