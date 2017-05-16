package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.models.ModelLoader;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.mappers.GenotypeFeatureMapper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.campagnelab.goby.util.Variant;
import org.campagnelab.goby.util.VariantMapHelper;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;
import java.util.Scanner;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */

public class DebugGenotype extends AbstractTool<DebugGenotypeArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    GenotypeFeatureMapper featureMapper;
    final public static boolean PRINT_INDEL_ERROR_CONTEXT = false;

    static private Logger LOG = LoggerFactory.getLogger(DebugGenotype.class);

    boolean queryFile = false;


    public static void main(String[] args) {

        DebugGenotype tool = new DebugGenotype();
        tool.parseArguments(args, "DebugGenotypeArguments", tool.createArguments());
        tool.execute();

    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {
        String genomePath = args().genomeFilename;

        try {
            if (args().featureMapperClassname != null) {
                    Class clazz = Class.forName(args().featureMapperClassname);
                    featureMapper = (GenotypeFeatureMapper) clazz.newInstance();
                    Properties sbiProp = new Properties();
                    sbiProp.load(new FileReader(args().inputFile+"p"));
                    featureMapper.configure(sbiProp);
            }
            System.err.println("Loading genome cache " + genomePath);
            genome = new RandomAccessSequenceCache();
            genome.load(genomePath, "min", "max");
            System.err.println("Done loading genome. ");
        } catch (ClassNotFoundException | IOException | IllegalAccessException | InstantiationException e) {
            System.err.println("Could not load genome cache");
            e.printStackTrace();
            System.exit(1);
        }
        VariantMapHelper varMap = null;
        RecordReader source = null;
        try {
            if (args().genotypeMap != null){
                varMap = new VariantMapHelper(args().genotypeMap);
            }
            source = new RecordReader(args().inputFile);
            queryFile = true;
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Unable to load true genotype map with filename " + args().genotypeMap, e);
        }

        Scanner sc = new Scanner(System.in);
        while (true) {
            System.out.println("Enter refId:pos, eg chr5:179535204 . enter \"exit\" to close. (1-indexed pos)");
            String input = sc.next();
            System.out.flush();
            if (input.equals("exit")){
                break;
            }
            String[] split = input.split(":");
            String chr = split[0];
            //goby pos
            int pos =  -1;
            try {
                pos = Integer.parseInt(split[1]) - 1;
            } catch (ArrayIndexOutOfBoundsException | NumberFormatException e){
                System.out.println("refId:pos badly formatted. try again or exit.");
                continue;
            }
            Integer pos2 = null;
            try {
                pos2 = Integer.parseInt(split[2]) - 1;
            } catch (ArrayIndexOutOfBoundsException | NumberFormatException e){
                continue;
            }

            //now find true genotype (map or genome)
            Variant var = null;
            if (varMap != null) {
                var = varMap.getVariant(chr,pos);
            } else {
                System.out.println("no varmap given");
            }
            String trueGenotype = "True Genotype ";
            if (var == null) {
                trueGenotype += " not in varmap:\n";
                int genomeTargetIndex = genome.getReferenceIndex(chr);
                if (genomeTargetIndex==-1 ) {
                    System.err.printf("Unable to locate reference sequence %s in genome.",chr);
                    break;
                }
                char referenceBaseChar = genome.get(genomeTargetIndex, pos);
                String referenceBase = Character.toString(referenceBaseChar);
                referenceBase = referenceBase.toUpperCase();
                trueGenotype += referenceBase + "|" + referenceBase;
            } else {
                trueGenotype += " found in map:\n";
                trueGenotype += var.toString();
            }
            System.out.println(trueGenotype);
            if (!queryFile){
                continue;
            }
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            System.out.println("scanning from total of " + source.numRecords() + " records.");
            for (BaseInformationRecords.BaseInformation rec :  source){
                if (pos2 != null){
                    if (rec.getPosition() >= pos && rec.getPosition() <= pos2 && rec.getReferenceId().equals(chr)){
                        System.out.println(rec.getReferenceId() + ":" + rec.getPosition() + ":" + rec.getSamples(0).getFormattedCounts());
                    }
                    if (rec.getPosition() == pos2){
                        break;
                    }
                    continue;
                }
                if (rec.getPosition() == pos && rec.getReferenceId().equals(chr)){
                    System.out.println ("\n SBI record: \n" + rec);
                    if (featureMapper!=null){
                        INDArray inputs = Nd4j.zeros(featureMapper.numberOfFeatures());
                        featureMapper.prepareToNormalize(rec,0);
                        featureMapper.mapFeatures(rec,inputs,0);
                        String[] featureVals = new String[6];
                        String[] featureNames = new String[6];
                        for (int i = 0; i < featureMapper.numberOfFeatures(); i++){
                            featureNames[i%6] = featureMapper.getFeatureName(i);
                            featureVals[i%6] = Float.toString(inputs.getFloat(i));
                            if (i % 6 == 5){
                                String names = String.format("%-30s %-30s %-30s %-30s %-30s", featureNames[0],
                                                                                                 featureNames[1],
                                                                                                 featureNames[2],
                                                                                                 featureNames[3],
                                                                                                 featureNames[4],
                                                                                                 featureNames[5]);
                                String vals = String.format("%-30s %-30s %-30s %-30s %-30s", featureVals[0],
                                                                                                 featureVals[1],
                                                                                                 featureVals[2],
                                                                                                 featureVals[3],
                                                                                                 featureVals[4],
                                                                                                 featureVals[5]);
                                System.out.println(names + "\n" + vals);
                                featureVals = new String[6];
                                featureNames = new String[6];
                            }
                        }

                    }

                    System.out.println(trueGenotype);
                    break;
                }
                recordLogger.lightUpdate();
            }
            recordLogger.done();

        }



    }


    @Override
    public DebugGenotypeArguments createArguments() {
        return new DebugGenotypeArguments();
    }


}
