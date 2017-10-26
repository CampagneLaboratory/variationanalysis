package org.campagnelab.dl.genotype.segments;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

/**
 * Util methods for {@link Segment}.
 */
public class SegmentUtil {

    /**
     * Determines if a genomic site has some reads that suggest an indel.
     *
     * @param baseInformation  the genomic site
     * @param candidateIndelThreshold a genomic site that has strictly more indel supporting reads than this threshold will be marked has candidateIndel.
     * @return
     */
    public static boolean hasCandidateIndel(BaseInformationRecords.BaseInformation baseInformation, int candidateIndelThreshold) {
        // determine if a genomic site has some reads that suggest an indel.
        for (BaseInformationRecords.SampleInfo sample : baseInformation.getSamplesList()) {
            for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
                if (count.getIsIndel()) {
                    if (count.getGenotypeCountForwardStrand() > candidateIndelThreshold || count.getGenotypeCountReverseStrand() > candidateIndelThreshold) {

                        return true;
                    }
                }
            }
        }
        return false;
    }

}
