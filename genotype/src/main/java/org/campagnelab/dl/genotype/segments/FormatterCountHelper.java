package org.campagnelab.dl.genotype.segments;

import it.unimi.dsi.lang.MutableString;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;



public class FormatterCountHelper {
    static public String format(BaseInformationRecords.SampleInfoOrBuilder sample){

        MutableString buffer = new MutableString();

        String from=null;
        for (BaseInformationRecords.CountInfo count : sample.getCountsList()) {
           if (from==null) {
               from=count.getFromSequence();
           }
            buffer.append(String.format(" %s=%d ", count.getToSequence(), count.getGenotypeCountForwardStrand() +
                    count.getGenotypeCountReverseStrand()));
        }

            buffer.append(" (from: " + from+")");
        return buffer.toString();
    }
}
