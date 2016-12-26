package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.logging.ProgressLogger;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Random;
import java.util.Set;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class AddTrueGenotypes extends AbstractTool<AddTrueGenotypesArguments> {


    int numVariantsAdded = 0;
    int numIndelsIgnored = 0;
    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();

    static private Logger LOG = LoggerFactory.getLogger(AddTrueGenotypes.class);

    public static void main(String[] args) {

        AddTrueGenotypes tool = new AddTrueGenotypes();
        tool.parseArguments(args, "AddTrueGenotypes", tool.createArguments());
        tool.execute();

    }


    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {


        //get reference genome
        String genomePath = args().genomeFilename;
        try {
            System.err.println("Loading genome cache " + genomePath);
            genome = new RandomAccessSequenceCache();
            genome.load(genomePath, "min", "max");
            System.err.println("Done loading genome. ");
        } catch (ClassNotFoundException e) {
            System.err.println("Could not load genome cache");
            e.printStackTrace();
            System.exit(1);
        } catch (IOException e) {
            System.err.println("Could not load genome cache");
            e.printStackTrace();
            System.exit(1);
        }

        int sampleIndex = args().sampleIndex;
        Random random = new XorShift1024StarRandom();
        try {
            Object2ObjectMap<String, Int2ObjectMap<String>> chMap = (Object2ObjectMap<String, Int2ObjectMap<String>>) BinIO.loadObject(args().genotypeMap);
            RecordReader source = new RecordReader(args().inputFile);
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);

            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            int recordsLabeled = 0;
            recordLogger.start();
            ObjectSet<String> distinctTrueGenotypes = new ObjectArraySet<>();
            for (BaseInformationRecords.BaseInformation rec : source) {
                boolean skip = false;


                BaseInformationRecords.BaseInformation.Builder buildRec = rec.toBuilder();
                int position = buildRec.getPosition();
                String chrom = buildRec.getReferenceId();
                Set<String> genotypeSet;
                String trueGenotype;
                int genomeTargetIndex = genome.getReferenceIndex(buildRec.getReferenceId());
                String referenceBase = Character.toString(genome.get(genomeTargetIndex, buildRec.getPosition()));
                boolean isVariant = false;
                boolean inMap = false;
                String genotypeFromMap = null;
                Int2ObjectMap<String> chromMap = chMap.get(chrom);

                if (chromMap != null) {
                    // The map contains Goby positions (zero-based).
                    genotypeFromMap = chromMap.get(position);
                    if (genotypeFromMap != null) {
                        inMap = true;
                    }
                }
                //indels should not be counted as variants: simply use refbase
                boolean isIndel = GenotypeHelper.isIndel(genotypeFromMap);
                if (isIndel) {
                    numIndelsIgnored++;
                }
                if (inMap && (!isIndel)) {
                    trueGenotype = genotypeFromMap;
                    if (!GenotypeHelper.isNoCall(genotypeFromMap)) {
                        isVariant = GenotypeHelper.isVariant(args().considerIndels /**/, genotypeFromMap, referenceBase);
                        if (isVariant) {
                            isVariant = true;
                            numVariantsAdded++;
                        }
                    }
                } else {
                    if (random.nextFloat() > args().referenceSamplingRate) {
                        skip = true;
                    }
                    // alignment and genome do not necessarily share the same space of reference indices. Convert:
                    trueGenotype = referenceBase + "|" + referenceBase;
                    referenceBase = referenceBase.toUpperCase();
                    genotypeSet = new ObjectArraySet<>();
                    genotypeSet.add(referenceBase);

                }

                if (!skip) {

                    distinctTrueGenotypes.add(trueGenotype);
                    // write the record.
                    buildRec.setTrueGenotype(trueGenotype.replace('|', '/').toUpperCase());
                    BaseInformationRecords.SampleInfo.Builder buildSample = buildRec.getSamples(sampleIndex).toBuilder();
                    for (int i = 0; i < buildSample.getCountsCount(); i++) {
                        BaseInformationRecords.CountInfo.Builder count = buildSample.getCounts(i).toBuilder();
                        boolean isCalled = GenotypeHelper.getAlleles(trueGenotype).contains(count.getToSequence());
                        count.setIsCalled(isCalled);
                        buildSample.setCounts(i, count);
                    }
                    buildSample.setIsVariant(isVariant);
                    buildRec.setSamples(sampleIndex, buildSample.build());
                    dest.appendEntry(buildRec.build());
                    recordsLabeled++;
                }
                recordLogger.lightUpdate();
            }
            recordLogger.done();
            dest.close();
            System.out.println("Found the following distinct true genotypes: " + distinctTrueGenotypes);
            System.out.println(numVariantsAdded + " number of variants in the sbi file.");
            System.out.println(numIndelsIgnored + " number of indels ignored in the file.");
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public AddTrueGenotypesArguments createArguments() {
        return new AddTrueGenotypesArguments();
    }


}
