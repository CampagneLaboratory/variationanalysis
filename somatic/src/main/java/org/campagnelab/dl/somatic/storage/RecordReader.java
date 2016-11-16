package org.campagnelab.dl.somatic.storage;

import org.apache.commons.io.IOUtils;
import org.campagnelab.dl.varanalysis.protobuf.BaseInformationRecords;
import org.campagnelab.goby.baseinfo.SequenceBaseInformationReader;
import org.campagnelab.goby.exception.GobyRuntimeException;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;
import java.util.Properties;
import java.util.Spliterator;
import java.util.function.Consumer;

/**
 * A reader for base information records stored in protobuf format.
 *
 * @author manuele simi
 */
public class RecordReader implements Closeable, RecordIterable, RecordReaderI<BaseInformationRecords.BaseInformation> {

    private SequenceBaseInformationReader reader;

    public RecordReader(String filepath) throws IOException {

        reader = new SequenceBaseInformationReader(filepath);
    }


    /**
     * Reads the next record, if available.
     *
     * @return the record
     * @throws IOException
     */
    public BaseInformationRecords.BaseInformation nextRecord() throws IOException {
        try {
            if (reader.hasNext()) return reader.next();
            else return null;
        } catch (GobyRuntimeException e) {return null;}
    }

    /**
     * Closes this stream and releases any system resources associated
     * with it. If the stream is already closed then invoking this
     * method has no effect.
     * <p>
     * <p> As noted in {@link AutoCloseable#close()}, cases where the
     * close may fail require careful attention. It is strongly advised
     * to relinquish the underlying resources and to internally
     * <em>mark</em> the {@code Closeable} as closed, prior to throwing
     * the {@code IOException}.
     *
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void close() throws IOException {
        IOUtils.closeQuietly(reader);
    }

    /**
     * Gets the number of records read so far.
     *
     * @return records loaded
     */
    public long getRecordsLoadedSoFar() {
        return reader.getRecordsLoadedSoFar();
    }

    /**
     * Gets the total number of records.
     *
     * @return total records
     */
    public long getTotalRecords() {
        return reader.getTotalRecords();
    }


    @Override
    public Iterator<BaseInformationRecords.BaseInformation> iterator() {
        return reader.iterator();
    }

    /**
     * Performs the given action for each element of the {@code Iterable}
     * until all elements have been processed or the action throws an
     * exception.  Unless otherwise specified by the implementing class,
     * actions are performed in the order of iteration (if an iteration order
     * is specified).  Exceptions thrown by the action are relayed to the
     * caller.
     *
     * @param action The action to be performed for each element
     * @throws NullPointerException if the specified action is null
     * @implSpec <p>The default implementation behaves as if:
     * <pre>{@code
     *     for (T t : this)
     *         action.accept(t);
     * }</pre>
     * @since 1.8
     */
    @Override
    public void forEach(Consumer<? super BaseInformationRecords.BaseInformation> action) {

    }

    /**
     * Creates a {@link Spliterator} over the elements described by this
     * {@code Iterable}.
     *
     * @return a {@code Spliterator} over the elements described by this
     * {@code Iterable}.
     * @implSpec The default implementation creates an
     * <em><a href="Spliterator.html#binding">early-binding</a></em>
     * spliterator from the iterable's {@code Iterator}.  The spliterator
     * inherits the <em>fail-fast</em> properties of the iterable's iterator.
     * @implNote The default implementation should usually be overridden.  The
     * spliterator returned by the default implementation has poor splitting
     * capabilities, is unsized, and does not report any spliterator
     * characteristics. Implementing classes can nearly always provide a
     * better implementation.
     * @since 1.8
     */
    @Override
    public Spliterator<BaseInformationRecords.BaseInformation> spliterator() {
        return null;
    }


    public Properties getProperties() {
        return reader.getProperties();
    }


    @Override
    public long numRecords() {
        return getTotalRecords();
    }
}
