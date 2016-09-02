package org.campagnelab.dl.varanalysis.intermediaries;

import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Iterator;

/**
 * Created by fac2003 on 7/19/16.
 * @author Fabien Campagne
 */
public class SimulationCharacteristics {


    ObjectArrayList<BaseInformationRecords.BaseInformation> records = new ObjectArrayList<>();

    public void observe(BaseInformationRecords.BaseInformation base) {
        records.add(base);
    }

    public int size() {
        return records.size();
    }

    public void clear() {
        records.clear();
    }

    public Iterator<BaseInformationRecords.BaseInformation> iterator() {
        return records.iterator();
    }

    /**
     * Observe characeteristics of bases on a complete batch.
     */
    public void batchIsComplete() {

    }
}
