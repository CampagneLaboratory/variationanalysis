package org.campagnelab.dl.genotype.helpers;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;

import java.io.IOException;
import java.util.Properties;
import java.util.Random;
import java.util.Set;

/**
 * Encapsulates the AddTrueGenotype logic.
 * Created by fac2003 on 12/27/16.
 */
public class AddTrueGenotypeHelper {
    private final RandomAccessSequenceCache genome;
    Random random = new XorShift1024StarRandom();
    Object2ObjectMap<String, Int2ObjectMap<String>> chMap;
    private int numIndelsIgnored;
    private int numVariantsAdded;
    private ObjectSet<String> distinctTrueGenotypes = new ObjectArraySet<>();
    private boolean considerIndels;
    private float referenceSamplingRate;
    private BaseInformationRecords.BaseInformation labeledEntry;
    private int sampleIndex;
    private int numRecords=0;
    private String mapFilename;

    /**
     * Create a helper with map, genome, etc.
     *
     * @param mapFilename
     * @param genome
     * @param sampleIndex
     * @param considerIndels
     * @param referenceSamplingRate
     */
    public AddTrueGenotypeHelper(String mapFilename, RandomAccessSequenceCache genome,
                                 int sampleIndex, boolean considerIndels, float referenceSamplingRate) {
        this.mapFilename=mapFilename;
        try {
            chMap = (Object2ObjectMap<String, Int2ObjectMap<String>>) BinIO.loadObject(mapFilename);
        } catch (IOException | ClassNotFoundException e) {

            throw new RuntimeException("Unable to load true genotype map with filename " + mapFilename, e);
        }
        this.genome = genome;
        this.considerIndels = considerIndels;
        this.referenceSamplingRate = referenceSamplingRate;
        this.sampleIndex = sampleIndex;
    }

    /**
     * Label a record with true genotype. Return true when the record should  be written to the output.
     * The result considers whether the site contains a true variant (always kept), or whether the site
     * has been sampled from reference matching sites (using referenceSamplingRate).
     *
     * @param record The .sbi record to annotate.
     * @return True if the record should be kept, i.e., written to the output, false otherwise.
     */
    public boolean addTrueGenotype(BaseInformationRecords.BaseInformation record) {

        numRecords++;

        boolean skip = false;
        BaseInformationRecords.BaseInformation.Builder buildRec = record.toBuilder();
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
                isVariant = GenotypeHelper.isVariant(considerIndels /**/, genotypeFromMap, referenceBase);
                if (isVariant) {
                    isVariant = true;
                    numVariantsAdded++;
                }
            }
        } else {
            if (random.nextFloat() > referenceSamplingRate) {
                skip = true;
            }
            // alignment and genome do not necessarily share the same space of reference indices. Convert:
            trueGenotype = referenceBase + "|" + referenceBase;
            referenceBase = referenceBase.toUpperCase();
            genotypeSet = new ObjectArraySet<>();
            genotypeSet.add(referenceBase);

        }
        final boolean keep = !skip;
        if (keep) {

            distinctTrueGenotypes.add(trueGenotype);
            // write the record.
            buildRec.setTrueGenotype(trueGenotype.replace('|', '/').toUpperCase());
            BaseInformationRecords.SampleInfo.Builder buildSample = buildRec.getSamples(sampleIndex).toBuilder();
            for (int i = 0; i < buildSample.getCountsCount(); i++) {
                BaseInformationRecords.CountInfo.Builder count = buildSample.getCounts(i).toBuilder();
                boolean isCalled = GenotypeHelper.genotypeHasAllele(trueGenotype, count.getToSequence());
                count.setIsCalled(isCalled);
                buildSample.setCounts(i, count);
            }
            buildSample.setIsVariant(isVariant);
            buildRec.setSamples(sampleIndex, buildSample.build());
            labeledEntry = buildRec.build();
        } else {
            labeledEntry = null;
        }
        return keep;
    }

    public BaseInformationRecords.BaseInformation labeledEntry() {
        return labeledEntry;
    }

    public int getNumIndelsIgnored() {
        return numIndelsIgnored;
    }

    public int getNumVariantsAdded() {
        return numVariantsAdded;
    }

    public Properties getStatProperties() {
        Properties result = new Properties();
        result.put("addTrueGenotypes.numIndelsIgnored", Integer.toString(numIndelsIgnored));
        result.put("addTrueGenotypes.numVariantsAdded", Integer.toString(numVariantsAdded));
        result.put("addTrueGenotypes.input.numRecords", Integer.toString(numRecords));
        result.put("addTrueGenotypes.referenceSamplingRate", Float.toString(referenceSamplingRate));
        result.put("addTrueGenotypes.considerIndels", Boolean.toString(considerIndels));
        result.put("addTrueGenotypes.mapFilename", mapFilename);
        return result;
    }
}
