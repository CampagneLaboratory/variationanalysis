package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;

import java.util.Iterator;


/**
 * Split result of the {@link SingleCandidateIndelSplitStrategy}
 */
public class SubSegment extends Segment {
    private final String candidateReferenceId;
    private final int candidateIndelPosition;
    private final long windowSize;

    protected SubSegment(final SingleCandidateIndelSplitStrategy.SizedBaseMap before,
                         final Segment parent, final BaseInformationRecords.BaseInformation indel,
                         long windowSize) {
        super(parent.fillInFeatures);
        this.candidateIndelPosition = indel.getPosition();
        this.candidateReferenceId = indel.getReferenceId();
        this.windowSize = windowSize;
        Iterator<BaseInformationRecords.BaseInformation> it = before.values().iterator();
        while (it.hasNext()) {
            this.add(it.next());
        }
        this.add(indel);
    }

    @Override
    public void add(BaseInformationRecords.BaseInformation base) {
        if (this.accept(base))
            super.add(base);
    }

    /**
     * Decides if the base belongs to this subsegment
     *
     * @param base
     * @return
     */
    private boolean accept(BaseInformationRecords.BaseInformation base) {
        if (base.getPosition() >= this.candidateIndelPosition)
            return base.getPosition() - this.candidateIndelPosition <= this.windowSize;
        else
            return this.candidateIndelPosition - base.getPosition() < this.windowSize;

    }

    /**
     * Checks if the subsegment has reached the window size
     *
     * @return
     */
    protected boolean isOpen() {
        return (this.getLastPosition() - this.candidateIndelPosition < this.windowSize);
    }


    public String getIndelReferenceId() {
        return this.candidateReferenceId;
    }

    public int getIndelPosition() {
        return this.candidateIndelPosition;
    }
}
