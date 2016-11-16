package org.campagnelab.dl.somatic.intermediaries;


import it.unimi.dsi.logging.ProgressLogger;
import org.campagnelab.dl.somatic.storage.DualReader;
import org.campagnelab.dl.somatic.storage.RecordWriter;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;

/**
 * Combines two parquet files, where one file was produced with Goby and filtering on, and the
 * other file with Goby with filtering turned off.
 * The files must be aligned by genomic position. The output contains the non-filtered counts in sample 0,
 * and the counts for filtered genotypes in sample 1. This file will be used to train a classifier to predict
 * genotypes that shoud be failed.
 *
 * @author Fabien Campagne
 */
public class FilteredNonFiltered {
    static int blockSize = 256 * 1024 * 1024;
    static int pageSize = 64 * 1024;
    static private Logger LOG = LoggerFactory.getLogger(FilteredNonFiltered.class);

    public static void main(String[] args) {

        //randomize
        FilteredNonFiltered tool = new FilteredNonFiltered();
        String inputFilenameA = args[0];
        String inputFilenameB = args[1];
        String outputFilename = args[2];
        new File(outputFilename).delete();
        tool.execute(inputFilenameA, inputFilenameB, outputFilename);
    }

    public void execute(String inputFilenameA, String inputFilenameB, String outputFilename) {
        try {
            DualReader reader = new DualReader(inputFilenameA, inputFilenameB);
            RecordWriter writer = new RecordWriter(outputFilename);
            ProgressLogger pg = new ProgressLogger(LOG);

            pg.expectedUpdates = reader.getTotalRecords();
            pg.start();
            while (reader.hasNext()) {
                BaseInformationRecords.BaseInformation nonFiltered = reader.first();
                BaseInformationRecords.BaseInformation filtered = reader.second();
                BaseInformationRecords.BaseInformation.Builder result = BaseInformationRecords.BaseInformation.newBuilder();
                result.mergeFrom(nonFiltered);
                BaseInformationRecords.SampleInfo.Builder sampleCopy;
                sampleCopy = BaseInformationRecords.SampleInfo.newBuilder(nonFiltered.getSamples(0));
                int genotypeIndex = 0;
                for (BaseInformationRecords.CountInfo.Builder countInfo :
                        sampleCopy.getCountsBuilderList()) {
                    int forwardFilteredCount = nonFiltered.getSamples(0).getCounts(genotypeIndex).getGenotypeCountForwardStrand() -
                            filtered.getSamples(0).getCounts(genotypeIndex).getGenotypeCountForwardStrand();
                    countInfo.setGenotypeCountForwardStrand(forwardFilteredCount);
                    int reverseFilteredCount = nonFiltered.getSamples(0).getCounts(genotypeIndex).getGenotypeCountReverseStrand() -
                            filtered.getSamples(0).getCounts(genotypeIndex).getGenotypeCountReverseStrand();
                    countInfo.setGenotypeCountReverseStrand(reverseFilteredCount);

                    genotypeIndex++;

                }
                result.addSamples(sampleCopy.build());
                // TextFormat.print(result, System.out);
                writer.writeRecord(result.build());
                pg.update();
            }

            pg.stop();
            reader.close();
            writer.close();

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
