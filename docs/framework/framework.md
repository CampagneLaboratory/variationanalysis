
Matcha provides an abstraction over [DeepLearning4J](https://deeplearning4j.org/) that aims to facilitate
the practical development and evaluation of deep learning models.

Since Matcha is a Java-based framework, it is best documented in Javadocs.
You can find the Javadocs at [javadoc.io](http://www.javadoc.io/doc/org.campagnelab.dl/framework)

This page provides pointers to help you understand the
 organizing principles of the framework. Use these pointers together
with the examples provided by the somatic module.

### Design
A deep learning project often requires all these elements:

 1. Convert data from some type of record to tensors, needed to train neural nets.
 The records must contain the features and labels the net will be trained with.
 Matcha is agnostic to the type of record, as long as you can represent it as an Object.
 We refer to it as ````RecordType````  below and in the code base.
 ````RecordType```` is a generic argument to most classes in the framework.
 2. Develop and optimize a neural network architecture to the specific problem
 3. Train models, often with early stopping on a validation set
 4. Evaluate model performance on additional test sets
 5. Use the model in production in other code

Matcha is designed to help you write small programs to help with each of these project steps.

### Vectorization
Before you can train a neural network, you need to organize training
data in tensors. Deeplearning4j provides the INDArray class to represent
 a tensor and your programs are responsible for populating tensors with data
  before starting training of a model.

#### Feature Mappers
 Matcha offers abstractions to help assemble tensors. It offers the ````FeatureMapper<RecordType>````
 abstraction whose responsibility is to accept data from a data record
  (often a line of text, or record of a structured data) and produce
  vectorized features.

#### Concatenating Mappers
Matcha also offers concat feature mappers, which accept two or more
FeatureMapper implementations and assembles the tensor that corresponds
to the concatenation of the features.  This gives the programmer the abilitity
to add or remove features without having to worry about maintaining
sound indexing in an INDArray (which can be an important source of bug
when using Deeplearning4j directly).

- [1D Concat](http://static.javadoc.io/org.campagnelab.dl/framework/1.1.1/org/campagnelab/dl/framework/mappers/ConcatFeatureMapper.html)


#### Label Mappers
The framework also offers label mappers, to produce vectorized labels for
training. The work similarly to feature mappers, but produce labels.
LabelMappers also have their concatenate mapper.

### Model Architecture
Matcha offers the concept of model architecture and makes it easy to swap
different model architectures to evaluate their impact on model performance.
Matcha support the DL4J ComputationGraph, allowing multiple inputs to
a neural net, and multiple outputs. Each input can be associated with its
own FeatureMapper and each output is associated with a LabelMapper.

A key difference between developing directly against DL4J and developing
with Matcha is that with Matcha, switching the model architecture can be
done by changing a command line argument when using the  train-model tool.

### DomainDescriptor
The [````DomainDescriptor<RecordType>````](http://static.javadoc.io/org.campagnelab.dl/framework/1.1.1/org/campagnelab/dl/framework/domains/DomainDescriptor.html) can be implemented for new domains to configure
which FeatureMappers, LabelMappers, ModelArchitecture and performance
metric should be used to train and evaluate a model.

The [````SomaticMutationDomainDescriptor````](https://github.com/CampagneLaboratory/variationanalysis/blob/master/somatic/src/main/java/org/campagnelab/dl/somatic/learning/domains/SomaticMutationDomainDescriptor.java), in the somatic module, provides a
complete example of implementation. Implementing the domain descriptor
is necessary before you can use Match to train models.

### Training a model
Once you have implemented FeatureMapper, LabelMapper and a DomainDescriptor,
you can subclass ````TrainModel<RecordType>````, the Matcha framework training tool.
This is easy to do as you can see in
[````TrainModelS````](https://github.com/CampagneLaboratory/variationanalysis/blob/master/somatic/src/main/java/org/campagnelab/dl/somatic/learning/TrainModelS.java)

Once you have completed this step, you get a training tool with useful
features to control the training of models for your project:
you can control the number of training records, adjust the learning rate
from the command line, adjust regularization or the amount of dropout.
Importantly, you get caching of mapped features and labels so that when
training on GPUs data starvation is greatly reduced and training can
proceed quickly.

### Predict on a test set
Implement a sub-class of  ````Predict<RecordType>```` and you are ready
to use the model to predict labels for a test set. See
the [PredictS](https://github.com/CampagneLaboratory/variationanalysis/blob/master/somatic/src/main/java/org/campagnelab/dl/somatic/tools/PredictS.java) implementation in the somatic module for an example.

### Show: a tool for error analysis
Implement a sub-class of  ````Show<RecordType>```` in order to visualize
records in custom formats. Show is designed to work with the output of
predict so that you can narrow down quickly on the most common errors
made by a model. Error analysis can help guide the design of better
features and model architectures.

### Keeping information organized
Matcha is designed to keep organized the information you produce in a project.

##### Models
Matcha writes trained models in a folder with all the information needed
to recreate the model from other Java code. Information is encoded in Java properties files
 to make it easy for both humans and machines to parse it.
 We even store the exact command line used to train the model.

##### Performance metric values
Performance meaured on the validation set is stored in the model-condition.txt
file at the end of training.  Using predict.sh will store performance in the predict-statistics.tsv


### Using models in other code

The somatic module is used in the Goby project to call somatic mutations.
Take a look at the following classes to see how this is done:

 - org.campagnelab.goby.modes.formats.SomaticVariationOutputFormat The class that uses the model in the Goby3 project.
 - org.campagnelab.goby.predictions.SomaticPredictor Interface in the Goby3 project.
 - org.campagnelab.dl.somatic.predictions.DLSomaticPredictor Implementation of the interface in the somatic module.
