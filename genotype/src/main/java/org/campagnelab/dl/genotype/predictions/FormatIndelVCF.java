package org.campagnelab.dl.genotype.predictions;

import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectArraySet;
import it.unimi.dsi.fastutil.objects.ObjectSet;
import org.campagnelab.dl.framework.domains.prediction.Prediction;
import org.campagnelab.dl.genotype.helpers.GenotypeHelper;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Collections;
import java.util.Set;
import java.util.SortedSet;

/**
 * Generate an indel compatible with VCF v4.1
 * Three steps:
 * 1. Set ref bases in to sequence set to the from sequence. This is because goby treats ref bases differently from VCF, doesn't extend them to match the indel's from sequence.
 *          IE: from: TGG to: T/T-G  -> from: TGG to: TGG/T-G
 * 2. trim all alleles to index of last dash any allele,
 *          IE: from: GTAC to: G--C/G-AC -> from: GTA to: G--/G-A
 * 3. delete dashes
 *          IE: from: GTA to: G--/G-A -> from: GTA to: G/GA
 * Created by rct66 on 2/7/16.
 */
public class FormatIndelVCF extends Prediction {

    public String fromVCF;
    public Set<String> toVCF;

    public FormatIndelVCF(String from, Set<String> to, char refBase){

        //step 1
        String refBaseStr = Character.toString(refBase);
        for (String s: to) {
            if (s.equals(refBaseStr)) {
                to.remove(s);
                to.add(from);
            }
        }


        //find newlen for step 2
        int newLen = -1;
        int maxLen = -1;
        for (String s : to){
            if (s.length() > maxLen){
                maxLen = s.length();
            }
        }
        maxLen = Math.max(maxLen,from.length());
        for (int i = 0; i < maxLen; i++){
            if (from.length() > i && from.charAt(i) == '-'){
                newLen = i+1;
            }
            for (String s : to){
                if (s.length() > i && s.charAt(i) == '-'){
                    newLen = i+1;
                }
            }
        }
        //apply step 2 and 3
        fromVCF = from.substring(0,Math.min(from.length(),newLen)).replace("-","");
        toVCF = new ObjectArraySet<>();
        for (String s : to){
            toVCF.add(s.substring(0,Math.min(s.length(),newLen)).replace("-",""));
        }




    }



}
