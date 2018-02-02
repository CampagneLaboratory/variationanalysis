package org.campagnelab.dl.framework.tools;

import org.deeplearning4j.nn.api.Model;
import org.deeplearning4j.nn.conf.NeuralNetConfiguration;
import org.deeplearning4j.nn.gradient.Gradient;
import org.deeplearning4j.optimize.api.ConvexOptimizer;
import org.deeplearning4j.optimize.api.IterationListener;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.primitives.Pair;

import java.util.Collection;
import java.util.Map;

public class PyTorchModel implements Model {
    private final String path;
    private final String label;

    public PyTorchModel(String path, String label) {
        this.path = path;
        this.label = label;
    }

    public String getPath() {
        return path;
    }

    public String getLabel() {
        return label;
    }

    @Override
    public void init() {

    }

    @Override
    public void setListeners(Collection<IterationListener> collection) {

    }

    @Override
    public void setListeners(IterationListener... iterationListeners) {

    }

    @Override
    public void addListeners(IterationListener... iterationListeners) {

    }

    @Override
    public void fit() {

    }

    @Override
    public void update(Gradient gradient) {

    }

    @Override
    public void update(INDArray indArray, String s) {

    }

    @Override
    public double score() {
        return 0;
    }

    @Override
    public void computeGradientAndScore() {

    }

    @Override
    public void accumulateScore(double v) {

    }

    @Override
    public INDArray params() {
        return null;
    }

    @Override
    public int numParams() {
        return 0;
    }

    @Override
    public int numParams(boolean b) {
        return 0;
    }

    @Override
    public void setParams(INDArray indArray) {

    }

    @Override
    public void setParamsViewArray(INDArray indArray) {

    }

    @Override
    public INDArray getGradientsViewArray() {
        return null;
    }

    @Override
    public void setBackpropGradientsViewArray(INDArray indArray) {

    }

    @Override
    public void applyLearningRateScoreDecay() {

    }

    @Override
    public void fit(INDArray indArray) {

    }

    @Override
    public void iterate(INDArray indArray) {

    }

    @Override
    public Gradient gradient() {
        return null;
    }

    @Override
    public Pair<Gradient, Double> gradientAndScore() {
        return null;
    }

    @Override
    public int batchSize() {
        return 0;
    }

    @Override
    public NeuralNetConfiguration conf() {
        return null;
    }

    @Override
    public void setConf(NeuralNetConfiguration neuralNetConfiguration) {

    }

    @Override
    public INDArray input() {
        return null;
    }

    @Override
    public void validateInput() {

    }

    @Override
    public ConvexOptimizer getOptimizer() {
        return null;
    }

    @Override
    public INDArray getParam(String s) {
        return null;
    }

    @Override
    public void initParams() {

    }

    @Override
    public Map<String, INDArray> paramTable() {
        return null;
    }

    @Override
    public Map<String, INDArray> paramTable(boolean b) {
        return null;
    }

    @Override
    public void setParamTable(Map<String, INDArray> map) {

    }

    @Override
    public void setParam(String s, INDArray indArray) {

    }

    @Override
    public void clear() {

    }
}
