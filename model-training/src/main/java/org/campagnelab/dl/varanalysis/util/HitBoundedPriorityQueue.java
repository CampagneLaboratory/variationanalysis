/*
 * Copyright (C) 2009-2010 Institute for Computational Biomedicine,
 *                    Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.campagnelab.dl.varanalysis.util;

import it.unimi.dsi.fastutil.objects.*;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;
import org.deeplearning4j.nn.multilayer.MultiLayerNetwork;
import org.nd4j.linalg.api.ndarray.INDArray;

import java.util.Collection;
import java.util.List;


/**
 * This class is an adaptation of MG4J's DocumentScoreBoundedSizeQueue.
 *
 * @author Fabien Campagne Date: Oct 21 2010
 */
public class HitBoundedPriorityQueue {

    /**
     * The underlying queue.
     */
    protected final ObjectSortedSet<ErrorRecord> queue;
    /**
     * The maximum number of documents to be ranked.
     */
    protected final int maxSize;


    /**
     * Creates a new empty bounded-size queue with a given capacity and natural
     * order as comparator. <p/> <P>Documents with equal scores will be compared
     * using their document index.
     *
     * @param capacity the initial capacity of this queue.
     */
    public HitBoundedPriorityQueue(final int capacity) {
        super();
        maxSize = capacity;
        queue = new ObjectAVLTreeSet<>(ErrorRecord.INCREASING_SCORE_COMPARATOR);

    }


    /**
     * Enqueues a transcript with given score and info.
     *
     * @return true if the document has been actually enqueued.
     */

    public synchronized boolean enqueue(float wrongness, INDArray features, INDArray label) {
        ErrorRecord dsi = new ErrorRecord(wrongness, features, label);
        if (queue.contains(dsi)) return false;
        if (queue.size() < maxSize) {
            queue.add(dsi);
            return true;
        } else {
            dsi = queue.first();

            if (wrongness > dsi.wrongness) {
                queue.remove(dsi);

                dsi.wrongness = wrongness;
                dsi.features = features;
                dsi.label = label;
                queue.add(dsi);
                return true;
            }
            return false;
        }
    }

    public boolean isEmpty() {
        return queue.isEmpty();
    }

    public int size() {
        return queue.size();
    }

    /**
     * Dequeues an error from the queue, returning an instance of {@link ErrorRecord}.
     * Documents are dequeued in inverse score order.
     *
     * @return the next {@link ErrorRecord}.
     */
    public final ErrorRecord dequeue() {
        ErrorRecord v = queue.first();
        queue.remove(v);
        return v;
    }

    public ErrorRecord first() {
        return queue.first();
    }

    public ErrorRecord last() {
        return queue.last();
    }

    public void clear() {
        queue.clear();
    }

    private XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom();

    /**
     * Return a random sample of exactly numElements.
     *
     * @param numElements The number of elements to return
     * @return a random sample drawn with replacement.
     */
    public List<ErrorRecord> getRandomSample(int numElements) {
        List<ErrorRecord> sample = new ObjectArrayList<>();
        assert !queue.isEmpty() : "queue cannot be empty to return a random sample";
        while (sample.size() < numElements) {
            for (ErrorRecord record : queue) {
                float value = random.nextFloat();
                if (value < 1f / queue.size()) {
                    sample.add(record);
                }
                if (sample.size() >= numElements) {
                    break;
                }
            }
        }
        return sample;
    }

    List<ErrorRecord> tmpList = new ObjectArrayList<>();

    public void updateWrongness(MultiLayerNetwork net) {
        tmpList.clear();

        for (ErrorRecord record : this.queue) {
            record.updateWrongness(record.features, net);
            tmpList.add(record);
        }
        //enqueue elements again to reorder them according to updated wrongness:
        queue.clear();
        queue.addAll(tmpList);

    }


    public float getMinWrongness() {
        return (float) queue.stream().mapToDouble((er -> er.wrongness)).min().getAsDouble();
    }

    public float getMaxWrongness() {
        return (float) queue.stream().mapToDouble((er -> er.wrongness)).max().getAsDouble();
    }

    public float getMeanWrongness() {
        return (float) queue.stream().mapToDouble((er -> er.wrongness)).sum() / queue.size();
    }
}
