package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.AddTrueGenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.campagnelab.goby.util.Variant;
import org.campagnelab.goby.util.VariantMapHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
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
        VariantMapHelper varMap = null;
        RecordReader source = null;
        try {
            if (args().genotypeMap != null){
                source = new RecordReader(args().inputFile);
                queryFile = true;
            }
            varMap = new VariantMapHelper(args().genotypeMap);
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
                break;
            }

            //now find true genotype (map or genome)
            Variant var = varMap.getVariant(chr,pos);
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
                if (rec.getPosition() == pos && rec.getReferenceId().equals(chr)){
                    System.out.println ("\n SBI record: \n" + rec);
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
