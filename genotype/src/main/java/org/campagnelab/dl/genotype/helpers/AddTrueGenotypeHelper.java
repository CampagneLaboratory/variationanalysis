package org.campagnelab.dl.genotype.helpers;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import it.unimi.dsi.util.XorShift1024StarRandom;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.predictions.AddTrueGenotypeHelperI;
import org.campagnelab.goby.reads.RandomAccessSequenceInterface;
import org.campagnelab.goby.util.Variant;
import org.campagnelab.goby.util.VariantMapHelper;
import org.campagnelab.goby.util.WarningCounter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.List;
import java.util.Properties;
import java.util.Random;
import java.util.Set;

/**
 * Encapsulates the AddTrueGenotype logic.
 * Created by fac2003 on 12/27/16.
 */
public class AddTrueGenotypeHelper implements AddTrueGenotypeHelperI {
    private final boolean SKIP_BAD_INDELS = true;

    private static final Logger LOG = LoggerFactory.getLogger(AddTrueGenotypeHelper.class);
    private RandomAccessSequenceInterface genome;
    Random random = new XorShift1024StarRandom();
    VariantMapHelper varMap;
    private int numIndelsIgnored;
    private int numIndelsAdded;
    private int numIndelsAddedAsRef;
    private int numSnpsAdded;
    private int numVariantsAdded;
    private int numHomozygousAdded;
    private int numHeterozygousAdded;
    private int numInMapAddedAsReference;
    private ObjectSet<String> distinctTrueGenotypes = new ObjectArraySet<>();
    private boolean considerIndels;
    private boolean indelsAsRef;
    private float referenceSamplingRate;
    private BaseInformationRecords.BaseInformation labeledEntry;
    private int sampleIndex;
    private int numRecords = 0;
    private int numWrongTrueCount = 0;
    private String mapFilename;
    private int recordsLabeled;
    static WarningCounter wrongNumGenosCalled = new WarningCounter(26);
    private List<BaseInformationRecords.BaseInformation> context;


    /**
     * Create a helper with map, genome, etc.
     *
     * @param mapFilename
     * @param genome
     * @param sampleIndex
     * @param considerIndels
     * @param referenceSamplingRate
     */
    public void configure(String mapFilename, RandomAccessSequenceInterface genome,
                          int sampleIndex, boolean considerIndels, boolean indelsAsRef, float referenceSamplingRate) {
        this.mapFilename = mapFilename;
        try {
            varMap = new VariantMapHelper(mapFilename);
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Unable to load true genotype map with filename " + mapFilename, e);
        }
        this.genome = genome;
        this.considerIndels = considerIndels;
        this.indelsAsRef = indelsAsRef;
        this.referenceSamplingRate = referenceSamplingRate;
        this.sampleIndex = sampleIndex;

    }

    public void configure(String mapFilename, RandomAccessSequenceInterface genome,
                          int sampleIndex, boolean considerIndels, float referenceSamplingRate) {
        this.mapFilename = mapFilename;
        try {
            varMap = new VariantMapHelper(mapFilename);
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Unable to load true genotype map with filename " + mapFilename, e);
        }
        this.genome = genome;
        this.considerIndels = considerIndels;
        this.indelsAsRef = true;
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
        int position = record.getPosition();
        String chrom = record.getReferenceId();
        int genomeTargetIndex = genome.getReferenceIndex(chrom);
        if (genomeTargetIndex==-1 ) {
            System.err.printf("Unable to locate reference sequence %s in genome.",chrom);
            System.exit(1);
        }
        char referenceBaseChar = genome.get(genomeTargetIndex, record.getPosition());
        String referenceBase = Character.toString(referenceBaseChar);
        return addTrueGenotype(willKeep(position, chrom, referenceBase), record);
    }

    // determine if a record will be kept
    public WillKeepI willKeep(int position, String referenceId, String referenceBase) {
        WillKeep willKepp = new WillKeep(position, referenceId, referenceBase).invoke();
        return willKepp;
    }

    //temporary addTrueGenotype that also receives context records for debugging
    public boolean addTrueGenotype(BaseInformationRecords.BaseInformation record, ObjectArrayList<BaseInformationRecords.BaseInformation> context) {
        this.context = context;
        return addTrueGenotype(record);
    }


    /**
     * Label a record with true genotype. Return true when the record should  be written to the output.
     * The result considers whether the site contains a true variant (always kept), or whether the site
     * has been sampled from reference matching sites (using referenceSamplingRate).
     *
     * @param willKeep Previously determined willKeep information.
     * @param record   The .sbi record to annotate.
     * @return True if the record should be kept, i.e., written to the output, false otherwise.
     */
    public boolean addTrueGenotype(WillKeepI willKeep, BaseInformationRecords.BaseInformation record) {

        numRecords++;
        // determine if the record should be kept:
        boolean keep = willKeep.isKeep();
        Set<Variant.FromTo> trueAlleles = willKeep.getTrueAlleles();
        String formattedTrueGenotype = GenotypeHelper.fromAlleles(GenotypeHelper.fromTosToAlleles(trueAlleles));
        boolean isVariant = willKeep.isVariant();
        BaseInformationRecords.BaseInformation.Builder buildRec = record.toBuilder();

        if (keep) {
            // We keep this record, so we label it:
            if (isVariant) {
                distinctTrueGenotypes.add(formattedTrueGenotype);
            }
            // write the record.
            buildRec.setTrueGenotype(formattedTrueGenotype);
            BaseInformationRecords.SampleInfo.Builder buildSample = buildRec.getSamples(sampleIndex).toBuilder();
            int trueAlleleCount = trueAlleles.size();
            StringBuilder matches = new StringBuilder("");
            for (int i = 0; i < buildSample.getCountsCount(); i++) {
                BaseInformationRecords.CountInfo.Builder count = buildSample.getCounts(i).toBuilder();
                String countFrom = count.getFromSequence();
                boolean isCalled = GenotypeHelper.genotypeHasAlleleOrIndel(trueAlleles,count.getToSequence(),countFrom);
                if (isCalled){
                    trueAlleleCount--;
                    matches.append(count.getFromSequence()).append(":").append(count.getToSequence()).append(", ");
                }
                count.setIsCalled(isCalled);
                buildSample.setCounts(i, count);
            }
            buildSample.setIsVariant(isVariant);
            buildRec.setSamples(sampleIndex, buildSample.build());
            labeledEntry = buildRec.build();
            recordsLabeled++;
            if (trueAlleleCount!=0){
//                StringBuffer sb = new StringBuffer();
//
//                sb.append(trueAlleleNum + " matching alleles expected\n");
//                sb.append("Too many or two few genotypes found: at Ref: " + record.getReferenceId() + " Pos: " + record.getPosition() + "\n" +
//                        "Ref:  " + trueFrom + " True:  " + trueGenotype + "\n");
//                for (BaseInformationRecords.CountInfo count: buildSample.getCountsList()){
//                    sb.append("from: " + count.getFromSequence() + " to: " + count.getToSequence() + " count: " + count.getGenotypeCountForwardStrand()+","+count.getGenotypeCountReverseStrand()+"\n");
//                }
//                sb.append("context: \n");
//                if (AddTrueGenotypes.PRINT_INDEL_ERROR_CONTEXT && context!=null){
//                    for (BaseInformationRecords.BaseInformation rec : context){
//                        sb.append("pos: " + rec.getPosition() + "\n");
//                        for (BaseInformationRecords.CountInfo count : rec.getSamples(0).getCountsList()){
//                            if (count.getIsIndel()){
//                                sb.append("from: " + count.getFromSequence() + " to: " + count.getToSequence() + " count: " + (count.getGenotypeCountForwardStrand() + count.getGenotypeCountReverseStrand()) + "\n");
//                            }
//                        }
//                    }
//                }
//
//                sb.append("Matches trimmedFrom:to : " + matches + "\n\n");
//                wrongNumGenosCalled.warn(LOG,sb.toString());
                if (SKIP_BAD_INDELS) {
                    keep = false;
                }
                numWrongTrueCount++;


            }
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

    public int getNumHomozygousAdded() {
        return numHomozygousAdded;
    }
    public int getNumInMapAddedAsReference() {
        return numInMapAddedAsReference;
    }
    public int getNumHeterozygousAdded() {
        return numHeterozygousAdded;
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

    public class WillKeep implements WillKeepI {

        private int position;
        private String chrom;

        private String referenceBase;
        private Set<Variant.FromTo> trueAlleles;
        private boolean isVariant;
        private boolean keep;
        private boolean isIndel;

        public WillKeep(int position, String chrom, String referenceBase) {

            this.position = position;
            this.chrom = chrom;
            this.referenceBase = referenceBase;
        }


        public boolean isVariant() {
            return isVariant;
        }

        public boolean isKeep() {
            return keep;
        }

        @Override
        public Set<Variant.FromTo> getTrueAlleles() {
            return trueAlleles;
        }

        public WillKeep invoke() {
            boolean isVariant = false;
            boolean skip = false;
            boolean inMap = false;
            boolean isIndel = false;
            // The map contains Goby positions (zero-based).
            Variant variant = varMap.getVariant(chrom,position);
            if (variant != null) {
                inMap = true;
            }
            if (inMap) {

                isIndel = variant.isIndel;
                trueAlleles = variant.trueAlleles;
                if (!GenotypeHelper.isNoCall(GenotypeHelper.fromAlleles(GenotypeHelper.fromTosToAlleles(variant.trueAlleles)))) {
                    isVariant = GenotypeHelper.isVariant(variant.trueAlleles);
                    if (isVariant) {
                        if (variant.isIndel) {
                            //indel is in map and is considered
                            numIndelsAdded++;
                        } else {
                            //snp is in map
                            numSnpsAdded++;
                        }
                        //we have a snp or indel
                        numVariantsAdded++;
                        if (trueAlleles.size() > 1) {
                            numHeterozygousAdded++;
                        } else {
                            numHomozygousAdded++;
                        }
                    } else {
                        //variant in map but not considered, to be added as ref
                        numInMapAddedAsReference++;
                    }
                }
                if (variant.isIndel && indelsAsRef && (!considerIndels)) {
                    //indel in map but added as ref
                    numIndelsAddedAsRef++;
                }
            }
            //handle case (whether or not is in map) where we want to use ref
            if (!isVariant) {
                if (random.nextFloat() > referenceSamplingRate) {
                    skip = true;
                }
                // alignment and genome do not necessarily share the same space of reference indices. Convert:
                referenceBase = referenceBase.toUpperCase();
                trueAlleles = new ObjectArraySet<Variant.FromTo>(1);
                trueAlleles.add(new Variant.FromTo(referenceBase,referenceBase));
            } else if (isVariant && variant.isIndel && (!indelsAsRef) && (!considerIndels)){
                numIndelsIgnored++;
                skip = true;
            }
            for (Variant.FromTo allele : trueAlleles){
                allele.makeUpperCase();
            }
            this.isVariant = isVariant;
            keep = !skip;
            return this;
        }
    }

    public void printStats() {

        int indelsSkipped = (SKIP_BAD_INDELS?numWrongTrueCount:0);

        System.out.println("Found the following distinct true genotypes (variants only): " + distinctTrueGenotypes);
        System.out.println(getNumVariantsAdded() - indelsSkipped + " number of variants in the sbi file.");
        System.out.println(getNumHeterozygousAdded() + " number of heterozygous added in the file.");
        System.out.println(getNumHomozygousAdded() + " number of homozygous added in the file.");
        System.out.println(getNumInMapAddedAsReference() + " number of map items added as reference instead.");
        System.out.println(recordsLabeled + " labeled records set for inclusion to file.");
        System.out.println(numIndelsAdded - indelsSkipped + " total number of indels added");
        System.out.println(numSnpsAdded + " total number of snps added");
        System.out.println(numIndelsAddedAsRef + " total number of indels added as ref");
        System.out.println(getNumIndelsIgnored() + " number of indels ignored instead of being added (as ref or indel)");
        System.out.println(numWrongTrueCount + " records found where too many or too few true genotypes set to called." + (SKIP_BAD_INDELS?" not included in sbi.":"included in sbi/"));

    }



}
