package org.campagnelab.dl.somatic.tools;


import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import org.campagnelab.dl.somatic.storage.RecordReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;

/**
 * The randomizer object iterates over a parquet file and randomizes the order of records in batches.
 * <p>
 * Created by rct66 on 5/18/16.
 *
 * @author rct66
 */
public class AddCalls extends AbstractTool<AddCallsArguments> {


    static private Logger LOG = LoggerFactory.getLogger(AddCalls.class);

    public static void main(String[] args) {

        AddCalls tool = new AddCalls();
        tool.parseArguments(args, "AddCalls", tool.createArguments());
        tool.execute();
    }


    @Override
    //only supports genotypes encoded with a bar (|) delimiter
    public void execute() {
        try {
            Object2ObjectMap<String,Int2ObjectMap<String>> chMap = (Object2ObjectMap<String,Int2ObjectMap<String>>)BinIO.loadObject(args().genotypeMap);
            RecordReader source = new RecordReader(args().inputFile);
            RecordWriter dest = new RecordWriter(args().inputFile+"_called");

            recordloop:
            for (BaseInformationRecords.BaseInformation rec : source) {
                //first, we skip examples where all reads match the reference
                for (BaseInformationRecords.CountInfo count : rec.getSamples(0).getCountsList()){
                    if (!count.getMatchesReference() && (count.getGenotypeCountForwardStrand() + count.getGenotypeCountReverseStrand()) != 0) {
                        break recordloop;
                    }
                }
                BaseInformationRecords.BaseInformation.Builder build = rec.toBuilder();
                int position = build.getPosition();
                String chrom = build.getReferenceId();
                String[] genotypes = chMap.get(chrom).get(position).split("|");
                for (BaseInformationRecords.CountInfo count : build.getSamples(0).getCountsList()){
                    BaseInformationRecords.CountInfo.Builder countBuilder = count.toBuilder();
                    boolean isCalled = (countBuilder.getToSequence().equals(genotypes[0])||countBuilder.getToSequence().equals(genotypes[1]));
                    countBuilder.setIsCalled(isCalled);
                }
                dest.writeRecord(build.build());
            }
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }



    @Override
    public AddCallsArguments createArguments() {
        return new AddCallsArguments();
    }



}
