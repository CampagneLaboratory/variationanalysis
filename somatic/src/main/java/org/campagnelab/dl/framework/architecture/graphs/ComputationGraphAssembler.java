package org.campagnelab.dl.framework.architecture.graphs;

import org.campagnelab.dl.framework.tools.TrainingArguments;
import org.campagnelab.dl.framework.domains.DomainDescriptor;
import org.campagnelab.dl.framework.models.ModelPropertiesHelper;
import org.deeplearning4j.nn.graph.ComputationGraph;
import org.nd4j.linalg.lossfunctions.LossFunctions;

/**
 * An interface for classes that assemble computation graphs with specific architectures.
 * A computation graph accepts a specific number of inputs, has a specific number of outputs,
 * and contains specific components (e.g., LSTM components, Dense layer component, etc.).
 * The graph may have more components than exposed via getComponentNames() if the number of nodes
 * of these components is not variable (fixed by the architecture).
 *
 * <p>Created by fac2003 on Nov 6 2016.</p>
 *
 * @author Fabien Campagne
 */
public interface ComputationGraphAssembler {
    /**
     * This method must be called before createComputationalGraph to give the assembler a chance to configure
     * the parameters described in the arguments (learning rate, dropout rate, etc.). The graph may have more
     * components than exposed via
     *
     * @param arguments Training arguments.
     */
    void setArguments(TrainingArguments arguments);

    /**
     * Create a new computational graph. Call this method after you configured the arguments, the number of floats in
     * each input and output of the computational graph.
     *
     * @return The fully configured computational graph, ready for training.
     */
    public ComputationGraph createComputationalGraph(DomainDescriptor domainDescriptor);

    /**
     * Set the dimensions of the specified input. You can describe 1-d inputs with only one dimension (the number of
     * floats expected in the input), 2-d inputs with two dimensions (e.g., useful for inputs to LSTM components,
     * first dimension is the number of features at each time point,  the second dimension is the number of time
     * points; also useful for images, where the first dimension can be the width and the second dimension the height).
     *
     * @param inputName name of the input for which we are setting dimensions.
     * @param dimension dimension.
     * @throws IllegalArgumentException when the outputName is not defined by this type of graph.
     */
    public void setNumInputs(String inputName, int... dimension);

    /**
     * Set the dimensions of the specified output. Similarly to inputs, you can describe 1-d, 2-d, etc dimensions.
     *
     * @param outputName name of the input for which we are setting dimensions.
     * @param dimension  dimension.
     * @throws IllegalArgumentException when the outputName is not defined by this type of graph.
     */
    public void setNumOutputs(String outputName, int... dimension);

    /**
     * Set the number of hidden nodes for a graph component.
     *
     * @param componentName  Name of the component.
     * @param numHiddenNodes Number of hidden nodes to configure for this component.
     * @throws IllegalArgumentException when the outputName is not defined by this type of graph.
     */
    public void setNumHiddenNodes(String componentName, int numHiddenNodes);

    /**
     * Returns the input names defined in this architecture.
     *
     * @return names of inputs.
     */
    String[] getInputNames();

    /**
     * Returns the output names defined in this architecture.
     *
     * @return names of outputs.
     */
    String[] getOutputNames();

    /**
     * Returns the component names defined in this architecture.
     *
     * @return names of components.
     */
    String[] getComponentNames();

    void setLossFunction(String outputName, LossFunctions.LossFunction lossFunction);

    /**
     * Save information about the specific architecture in the model properties.
     * @param helper
     */
    void saveProperties(ModelPropertiesHelper helper);
}
