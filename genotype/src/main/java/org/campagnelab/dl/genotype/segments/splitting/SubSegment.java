package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Iterator;

/**
 * Created by mas2182 on 10/25/17.
 */
public class SubSegment {
    private final Segment segment;
    private final int candidateIndelPosition;

    protected SubSegment(SingleCandidateIndelSplitStrategy.SizedBaseMap before, Segment parent, BaseInformationRecords.BaseInformation indel) {
           // synchronized (before) {
                Iterator<BaseInformationRecords.BaseInformation> it = before.values().iterator();
                segment = new Segment(parent.fillInFeatures, it.next());
                while (it.hasNext()) {
                    BaseInformationRecords.BaseInformation base = it.next();
                    segment.add(base);
                }
                segment.add(indel);

            //}
            this.candidateIndelPosition = indel.getPosition();
    }

    protected void addToAfterList(BaseInformationRecords.BaseInformation base) {
        
    }

    protected boolean isClosed() {
        return false;
    }

    protected Segment asSegment() {
        return this.segment;
    }
}
