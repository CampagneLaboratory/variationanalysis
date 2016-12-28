package org.campagnelab.dl.framework.mappers.processing;

/**
 * Helper class to facilitate setting bounds on masked elements in processing mappers
 *
 * Created by joshuacohen on 12/13/16.
 */
class Bounds {
    private int start;
    private int end;
    private int shiftedSize;

    Bounds() {
        shiftedSize = 0;
    }

    void setStart(int start) {
        this.start = start;
    }

    void setEnd(int end) {
        this.end = end;
    }

    void setShiftedSize(Bounds otherBounds) {
        shiftedSize = otherBounds.shiftedSize + otherBounds.size();
        start -= shiftedSize;
        end -= shiftedSize;
    }

    boolean contains(int index) {
        return start <= index;
    }

    int size() {
        return end - start;
    }
}
