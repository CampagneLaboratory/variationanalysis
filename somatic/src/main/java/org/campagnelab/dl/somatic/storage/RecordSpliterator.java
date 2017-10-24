package org.campagnelab.dl.somatic.storage;

import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.compression.ChunkCodec;
import org.campagnelab.goby.compression.MessageChunksWriter;

import java.io.File;
import java.io.IOException;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Consumer;

/**
 * Spliterator for {@link BaseInformationRecords.BaseInformation}
 *
 * @author manuele
 */
public class RecordSpliterator implements Spliterator<BaseInformationRecords.BaseInformation> {

    private final RecordReader reader;
    private final long startIndex; // the index of the first byte to read
    private final String sourceSBI;
    private long endIndex; // the index of the last byte to read
    private final long maxLength; //length in bytes of the input SBI
    private final long minLength; // minimun number of bytes assigned to a single spliterator

    public RecordSpliterator(String sourceSBI, long startIndex, long endIndex, long maxLength) {
        this.sourceSBI = sourceSBI;
        try {
            this.reader = new RecordReader(sourceSBI,startIndex,endIndex);
        } catch (IOException e) {
            throw new RuntimeException("unable to open the source sbi " + sourceSBI);
        }
        this.maxLength = maxLength;
        this.minLength = (maxLength / 10);
        this.startIndex = startIndex;
        this.endIndex = endIndex;
    }

    /**
     * If a remaining element exists, performs the given action on it,
     * returning {@code true}; else returns {@code false}.  If this
     * Spliterator is {@link #ORDERED} the action is performed on the
     * next element in encounter order.  Exceptions thrown by the
     * action are relayed to the caller.
     *
     * @param action The action
     * @return {@code false} if no remaining elements existed
     * upon entry to this method, else {@code true}.
     * @throws NullPointerException if the specified action is null
     */
    @Override
    public boolean tryAdvance(Consumer<? super BaseInformationRecords.BaseInformation> action) {
        if (action == null)
            throw new NullPointerException();
        BaseInformationRecords.BaseInformation record = null;
        try {
            record = reader.nextRecord();
        } catch (IOException e) {
            return false;
        }
        if (record != null) {
            action.accept(record);
            return true;
        }
        return false;

    }

    /**
     * If this spliterator can be partitioned, returns a Spliterator
     * covering elements, that will, upon return from this method, not
     * be covered by this Spliterator.
     * <p>
     * <p>If this Spliterator is {@link #ORDERED}, the returned Spliterator
     * must cover a strict prefix of the elements.
     * <p>
     * <p>Unless this Spliterator covers an infinite number of elements,
     * repeated calls to {@code trySplit()} must eventually return {@code null}.
     * Upon non-null return:
     * <ul>
     * <li>the value reported for {@code estimateSize()} before splitting,
     * must, after splitting, be greater than or equal to {@code estimateSize()}
     * for this and the returned Spliterator; and</li>
     * <li>if this Spliterator is {@code SUBSIZED}, then {@code estimateSize()}
     * for this spliterator before splitting must be equal to the sum of
     * {@code estimateSize()} for this and the returned Spliterator after
     * splitting.</li>
     * </ul>
     * <p>
     * <p>This method may return {@code null} for any reason,
     * including emptiness, inability to split after traversal has
     * commenced, data structure constraints, and efficiency
     * considerations.
     *
     * @return a {@code Spliterator} covering some portion of the
     * elements, or {@code null} if this spliterator cannot be split
     * @apiNote An ideal {@code trySplit} method efficiently (without
     * traversal) divides its elements exactly in half, allowing
     * balanced parallel computation.  Many departures from this ideal
     * remain highly effective; for example, only approximately
     * splitting an approximately balanced tree, or for a tree in
     * which leaf nodes may contain either one or two elements,
     * failing to further split these nodes.  However, large
     * deviations in balance and/or overly inefficient {@code
     * trySplit} mechanics typically result in poor parallel
     * performance.
     */
    @Override
    public Spliterator<BaseInformationRecords.BaseInformation> trySplit() {
        if ( ((endIndex - startIndex < 0) || (endIndex - startIndex < 100000000)))
            return null;
        long splitPosition = ((endIndex - startIndex) / 2) + this.startIndex; //we try to split halfway
        long readPosition = this.reader.getCurrentReadPosition();
        if (splitPosition < readPosition ) {
            //if we passed the current position of the reader, we don't further split
            return null;
        }
        Spliterator<BaseInformationRecords.BaseInformation> newSp = new RecordSpliterator(sourceSBI,splitPosition,endIndex,maxLength);
        this.endIndex = splitPosition - 1;
        this.reader.readUpTo(endIndex);
        return newSp;
    }

    private void recalculateEndIndex() {
        this.endIndex = startIndex + Math.round(this.endIndex / 2);
    }

    /**
     * Returns an estimate of the number of elements that would be
     * encountered by a {@link #forEachRemaining} traversal, or returns {@link
     * Long#MAX_VALUE} if infinite, unknown, or too expensive to compute.
     * <p>
     * <p>If this Spliterator is {@link #SIZED} and has not yet been partially
     * traversed or split, or this Spliterator is {@link #SUBSIZED} and has
     * not yet been partially traversed, this estimate must be an accurate
     * count of elements that would be encountered by a complete traversal.
     * Otherwise, this estimate may be arbitrarily inaccurate, but must decrease
     * as specified across invocations of {@link #trySplit}.
     *
     * @return the estimated size, or {@code Long.MAX_VALUE} if infinite,
     * unknown, or too expensive to compute.
     * @apiNote Even an inexact estimate is often useful and inexpensive to compute.
     * For example, a sub-spliterator of an approximately balanced binary tree
     * may return a value that estimates the number of elements to be half of
     * that of its parent; if the root Spliterator does not maintain an
     * accurate count, it could estimate size to be the power of two
     * corresponding to its maximum depth.
     */
    @Override
    public long estimateSize() {
        //how many bytes we would read
        long est = (this.endIndex - this.startIndex) * this.reader.getTotalRecords();
        return est>0 ? est : 0L;
    }

    /**
     * Returns a set of characteristics of this Spliterator and its
     * elements. The result is represented as ORed values from {@link
     * #ORDERED}, {@link #DISTINCT}, {@link #SORTED}, {@link #SIZED},
     * {@link #NONNULL}, {@link #IMMUTABLE}, {@link #CONCURRENT},
     * {@link #SUBSIZED}.  Repeated calls to {@code characteristics()} on
     * a given spliterator, prior to or in-between calls to {@code trySplit},
     * should always return the same result.
     * <p>
     * <p>If a Spliterator reports an inconsistent set of
     * characteristics (either those returned from a single invocation
     * or across multiple invocations), no guarantees can be made
     * about any computation using this Spliterator.
     *
     * @return a representation of characteristics
     * @apiNote The characteristics of a given spliterator before splitting
     * may differ from the characteristics after splitting.  For specific
     * examples see the characteristic values {@link #SIZED}, {@link #SUBSIZED}
     * and {@link #CONCURRENT}.
     */
    @Override
    public int characteristics() {
        return ORDERED | IMMUTABLE | NONNULL;
    }
}
