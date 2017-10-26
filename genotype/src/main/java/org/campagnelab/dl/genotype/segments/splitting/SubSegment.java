package org.campagnelab.dl.genotype.segments.splitting;

import org.campagnelab.dl.genotype.segments.Segment;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.dl.varanalysis.protobuf.SegmentInformationRecords;

import java.util.Iterator;
import java.util.Objects;
import java.util.function.Function;

/**
 * Split result of the {@link SingleCandidateIndelSplitStrategy}
 */
public class SubSegment {
    private Segment delegate = null;
    private final int candidateIndelPosition;
    private final long windowSize;

    protected SubSegment(final SingleCandidateIndelSplitStrategy.SizedBaseMap before,
                         final Segment parent, final BaseInformationRecords.BaseInformation indel,
                         long windowSize) {
        this.candidateIndelPosition = indel.getPosition();
        this.windowSize = windowSize;
        Iterator<BaseInformationRecords.BaseInformation> it = before.values().iterator();
        while (it.hasNext()) {
            this.addToOrCreateSegment(parent.fillInFeatures, it.next());
        }
        this.addToOrCreateSegment(parent.fillInFeatures, indel);
        Objects.nonNull(delegate);
    }

    private void addToOrCreateSegment(Function<BaseInformationRecords.BaseInformation, SegmentInformationRecords.Base.Builder> fillInFeatures,
                                      BaseInformationRecords.BaseInformation base) {
        if (this.accept(base)) {
            if (delegate == null)
                delegate = new Segment(fillInFeatures, base);
            else
                this.delegate.add(base);
        }
    }

    protected void add(BaseInformationRecords.BaseInformation base) {
        if (this.accept(base))
            this.delegate.add(base);
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
        return (this.delegate.getLastPosition() - this.candidateIndelPosition < this.windowSize);
    }

    protected Segment asSegment() {
        return this.delegate;
    }

    public int getIndelPosition() {
        return this.candidateIndelPosition;
    }

    public int getFirstPosition() {
        return this.delegate.getFirstPosition();
    }

    public int getLastPosition() {
        return this.delegate.getLastPosition();
    }
}
