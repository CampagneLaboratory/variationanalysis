package org.campagnelab.dl.genotype.tools;


import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.logging.ProgressLogger;
import org.apache.commons.io.FilenameUtils;
import org.campagnelab.dl.framework.tools.arguments.AbstractTool;
import org.campagnelab.dl.genotype.helpers.AddTrueGenotypeHelper;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationWriter;
import org.campagnelab.goby.reads.RandomAccessSequenceCache;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * The addcalls object uses a map to create a new protobuf file with genotype calls.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class AddTrueGenotypesMultiThreaded extends AbstractTool<AddTrueGenotypesArguments> implements Runnable{


    RandomAccessSequenceCache genome;
    final static int THREAD_COUNT = 8;
    final public static boolean PRINT_INDEL_ERROR_CONTEXT = false;

    static private Logger LOG = LoggerFactory.getLogger(AddTrueGenotypesMultiThreaded.class);

    public static void main(String[] args) {

        AddTrueGenotypesMultiThreaded tool = new AddTrueGenotypesMultiThreaded();
        tool.parseArguments(args, "AddTrueGenotypes", tool.createArguments());
        tool.prepare();
        for (int i = 0; i < THREAD_COUNT; i++ ){
            new Thread(tool).start();
        }
        tool.finish();

    }

    public void prepare(){
        try {
            source = new RecordReader(args().inputFile);
            iter =  source.iterator();
            recordLogger = new ProgressLogger(LOG);
            recordLogger.expectedUpdates = source.numRecords();
            System.out.println(source.numRecords() + " records to label");
            recordsLabeled = 0;
            recordLogger.start();
            ObjectArrayList<BaseInformationRecords.BaseInformation> recContext = new ObjectArrayList<>(1000);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void run() {
        AddTrueGenotypeHelper addTrueGenotypeHelper  = new AddTrueGenotypeHelper();
        SequenceBaseInformationWriter dest = null;
        genome = new RandomAccessSequenceCache();
        try {
            genome.load(args().genomeFilename, "min", "max");
            dest = new SequenceBaseInformationWriter(FilenameUtils.getFullPath(args().outputFilename) + FilenameUtils.getBaseName(args().outputFilename) + Thread.currentThread().getId() + ".sbi");
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        synchronized (this) {
            writers.add(dest);
        }
        addTrueGenotypeHelper.configure(
                args().genotypeMap,
                genome,
                args().sampleIndex,
                args().considerIndels,
                args().indelsAsRef,
                args().referenceSamplingRate);
        try {
            BaseInformationRecords.BaseInformation rec;
            boolean keep;
            while (true){
                synchronized(this) {
                    if (!iter.hasNext()){
                        break;
                    }
                    rec = iter.next();
                }
                keep = addTrueGenotypeHelper.addTrueGenotype(rec);
                synchronized(this) {
                    if (keep) {
                        dest.appendEntry(addTrueGenotypeHelper.labeledEntry());
                    }
                    recordsLabeled++;
                    recordLogger.lightUpdate();
                }



            }
            dest.setCustomProperties(addTrueGenotypeHelper.getStatProperties());
            dest.close();
            synchronized (this) {
                addTrueGenotypeHelper.printStats();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    RecordReader source;
    ProgressLogger recordLogger;

    int recordsLabeled;
    Iterator<BaseInformationRecords.BaseInformation> iter;
    List<SequenceBaseInformationWriter> writers = new ArrayList<SequenceBaseInformationWriter>();




    public void finish(){
        recordLogger.done();
        for(SequenceBaseInformationWriter writer : writers){
            System.out.println(writer.toString());
        }


    }

    @Override
    public AddTrueGenotypesArguments createArguments() {
        return new AddTrueGenotypesArguments();
    }

    @Override
    public void execute() {
        //use run instead
    }

}
