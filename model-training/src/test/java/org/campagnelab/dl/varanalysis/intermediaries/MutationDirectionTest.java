package org.campagnelab.dl.varanalysis.intermediaries;

import com.google.protobuf.TextFormat;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.junit.Test;

import java.io.StringWriter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * Test the mutator on some specific examples.
 * Created by fac2003 on 5/27/16.
 */
public class MutationDirectionTest {

    @Test
    public void mutateTest() throws Exception {
        FirstSimulationStrategy strategy = new FirstSimulationStrategy(1,1,0.1,1);
        final int mostCount = 100;
        final int secondMostCount = 11;
        //homozygous first:
        for (int numGenos = 5; numGenos < 8; numGenos++){
            for (int maxIndex = 0; maxIndex < numGenos; maxIndex++){
                if (maxIndex == 4){
                    break;
                }
                int[] counts = new int[numGenos];
                counts[maxIndex] = mostCount;
                for (int trial = 0; trial < 10; trial++){
                    FirstSimulationStrategy.mutationDirection dir = strategy.dirFromCounts(counts);
                    assert(dir.newBase !=dir.oldBase);
                    assert(dir.newBase !=maxIndex);
                    assert(dir.newBase != 4);
                    assert(dir.newBase < numGenos);
                }
            }
        }
        //now heterozygous
        for (int numGenos = 5; numGenos < 8; numGenos++){
            for (int maxIndex = 0; maxIndex < numGenos; maxIndex++){
                if (maxIndex == 4){
                    break;
                }
                for (int secondMostIndex = 0; secondMostIndex < numGenos; secondMostIndex++){
                    if (secondMostIndex == maxIndex || secondMostIndex == 4){
                        break;
                    }
                    int[] counts = new int[numGenos];
                    counts[maxIndex] = mostCount;
                    counts[secondMostIndex] = secondMostCount;
                    for (int trial = 0; trial < 10; trial++){
                        FirstSimulationStrategy.mutationDirection dir = strategy.dirFromCounts(counts);
                        assert(dir.newBase != dir.oldBase);
                        assert(dir.newBase != maxIndex);
                        assert(dir.newBase != 4);
                        assert(dir.newBase != secondMostIndex);
                        assert(dir.newBase < numGenos);
                    }
                }

            }
        }

    }


}