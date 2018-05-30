package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.bed.BEDRecords;
import org.campagnelab.dl.framework.bed.BedLoader;
import org.campagnelab.dl.framework.bed.FullOverlapper;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.AddTrueGenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileReader;
import java.io.IOException;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class AddTrueGenotypes extends AbstractTool<AddTrueGenotypesArguments> {


    RandomAccessSequenceCache genome = new RandomAccessSequenceCache();
    final public static boolean PRINT_INDEL_ERROR_CONTEXT = false;

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
        BEDRecords confidenceRegions = new FullOverlapper();
        if (args().confidenceRegionsFilename != null) {
            System.out.println("Loading confidence regions from " + args().confidenceRegionsFilename);
            try {
                confidenceRegions = BedLoader.loadBedFile(new FileReader(args().confidenceRegionsFilename));
            } catch (IOException e) {
                System.err.println("Unable to load BED file with path " + args().confidenceRegionsFilename);
                System.exit(1);
            }
            System.out.printf("Done loading %d confidence regions. Output will be restricted to sites in these regions.",
                    confidenceRegions.numRecords());
        }
        try {
            RecordReader source = new RecordReader(args().inputFile);
            SequenceBaseInformationWriter dest = new SequenceBaseInformationWriter(args().outputFilename);
            AddTrueGenotypeHelper addTrueGenotypeHelper = new AddTrueGenotypeHelper();
            addTrueGenotypeHelper.configure(
                    args().genotypeMap,
                    genome,
                    args().sampleIndex,
                    args().considerIndels,
                    args().indelsAsRef,
                    args().referenceSamplingRate);
            ProgressLogger recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            int recordsLabeled = 0;
            recordLogger.start();
            int doesNotOverlapConfidenceRegionsCount = 0;
            ObjectArrayList<BaseInformationRecords.BaseInformation> recContext = new ObjectArrayList<>(1000);
            for (BaseInformationRecords.BaseInformation rec : source) {
                boolean keep = false;
                final boolean overlapsConfidenceRegions = confidenceRegions.overlaps(rec.getReferenceId(), rec.getPosition(), rec.getPosition() + 1);
                doesNotOverlapConfidenceRegionsCount += (overlapsConfidenceRegions == false) ? 1 : 0;
                if (PRINT_INDEL_ERROR_CONTEXT) {
                    recContext.add(rec);
                    if (recContext.size() < 50) {
                        continue;
                    }
                    keep = addTrueGenotypeHelper.addTrueGenotype(recContext.get(recContext.size() / 2), recContext);
                    // filter by confidence region if provided (if not we use the fullOverlapper that does not change the keep flag):
                    keep &= overlapsConfidenceRegions;
                    if (keep) {
                        dest.appendEntry(addTrueGenotypeHelper.labeledEntry());
                    }

                    recContext.remove(0);
                } else {
                    keep = addTrueGenotypeHelper.addTrueGenotype(rec);
                    // filter by confidence region if provided (if not we use the fullOverlapper that does not change the keep flag):
                    keep &= overlapsConfidenceRegions;

                    if (keep) {
                        dest.appendEntry(addTrueGenotypeHelper.labeledEntry());
                    }
                }
                recordsLabeled++;
                recordLogger.lightUpdate();


            }
            recordLogger.done();
            dest.setCustomProperties(addTrueGenotypeHelper.getStatProperties());
            dest.close();
            addTrueGenotypeHelper.printStats();
            if (args().confidenceRegionsFilename != null) {
                System.out.printf("Removed %d sites outside of confidence regions.",doesNotOverlapConfidenceRegionsCount);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }


    @Override
    public AddTrueGenotypesArguments createArguments() {
        return new AddTrueGenotypesArguments();
    }


}
