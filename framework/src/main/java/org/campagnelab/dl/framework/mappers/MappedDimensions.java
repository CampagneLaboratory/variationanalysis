package org.campagnelab.dl.framework.mappers;

import java.util.Arrays;

/**
 * Information about the dimensions of a tensor.
 *
 * @author Fabien Campagne
 */
public class MappedDimensions {
    public int[] dimensions;

    /**
     * Create a new mapped dimension instance. The dimension values are the number of elements that can fit in each
     * dimension.
     *
     * @param dimension One or more dimensions of a tensor. 1 dimension when 1d, 2 dimensions when 2d, etc.
     */
    public MappedDimensions(int... dimension) {
        this.dimensions = dimension;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(dimensions);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof MappedDimensions) {
            MappedDimensions other = (MappedDimensions) obj;
            return Arrays.equals(other.dimensions, dimensions);
        }
        return false;
    }

    /**
     * Compare dimension at dimIdx to same dimension in other MappedDimensions
     *
     * @param other the other MappedDimensions instance to compare against
     * @param dimIdx the index of the dimension (1 for first dimension, 2 for 2nd, etc)
     * @return true if the dimension at dimIdx equals the same dimension in other
     */
    public boolean equalsDimension(MappedDimensions other, int dimIdx) {
        if (dimIdx < 1 || dimIdx > dimensions.length) {
            throw new IllegalArgumentException("Invalid dimension provided");
        }
        return other.dimensions[dimIdx - 1] == dimensions[dimIdx - 1];
    }

    /**
     * Return the number of elements in a tensor with these dimensions (product of the dimensions).
     * A 1-d tensor with the first dimension=30 will
     * return 30. A 2-d tensor with dimensions 10,2 will return 20, and so on.
     *
     * @return Number of elements.
     */
    public int numElements() {
        int num = 1;
        for (int dim : dimensions) {
            num *= dim;
        }
        return num;
    }

    /**
     * Return the number of elements in a tensor with these dimensions, at the dimension specified
     *
     * @param dimIdx the index of the dimension (1 for first dimension, 2 for 2nd, etc)
     * @return Number of elements at dimIdx
     */
    public int numElements(int dimIdx) {
        if (dimIdx < 1 || dimIdx > dimensions.length) {
            throw new IllegalArgumentException("Invalid dimension provided");
        }
        return dimensions[dimIdx - 1];
    }

    /**
     * The number of dimensions in the tensor.
     *
     * @return number of dimensions, 1 for 1d, 2 for 2d, etc.
     */
    public int numDimensions() {
        return dimensions.length;
    }
}
