//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef LIBTORCH_ENABLED

#include "LibtorchNeuralNetBase.h"
#include "DataIO.h"
#include "MultiMooseEnum.h"

namespace Moose
{

// A class that describes a simple feed-forward neural net.
class LibtorchArtificialNeuralNet : public LibtorchNeuralNetBase
{
public:
  /**
   * Construct using input parameters
   * @param name Name of the neural network
   * @param num_inputs The number of input neurons/parameters
   * @param num_neurons_per_layer Number of neurons per hidden layer
   * @param num_outputs The number of output neurons
   */
  LibtorchArtificialNeuralNet(const std::string name,
                              const unsigned int num_inputs,
                              const unsigned int num_outputs,
                              const std::vector<unsigned int> & num_neurons_per_layer,
                              const std::vector<std::string> & activation_function = {"relu"});

  /**
   * Add layers to the neural network
   * @param layer_name The name of the layer to be added
   * @param parameters A map of parameter names and the corresponding values which
   *                   describe the neural net layer architecture
   */
  virtual void addLayer(const std::string & layer_name,
                        const std::unordered_map<std::string, unsigned int> & parameters) override;

  /**
   * Overriding the forward substitution function for the neural network, unfortunately
   * this cannot be const since it creates a graph in the background
   * @param x Input tensor for the evaluation
   */
  torch::Tensor forward(torch::Tensor x);

  /// Return the name of the neural network
  const std::string & name() const { return _name; }
  /// Return the number of neurons on the input layer
  unsigned int numInputs() const { return _num_inputs; }
  /// Return the number of neurons on the output layer
  unsigned int numOutputs() const { return _num_outputs; }
  /// Return the number of hidden layers
  unsigned int numHiddenLayers() const { return _num_hidden_layers; }
  /// Return the hidden layer architecture
  const std::vector<unsigned int> & numNeuronsPerLayer() const { return _num_neurons_per_layer; }
  /// Return the multi enum containing the activation functions
  const MultiMooseEnum & activationFunctions() const { return _activation_function; }

  /// Construct the neural network
  void constructNeuralNetwork();

protected:
  /// Name of the neural network
  const std::string _name;
  /// Submodules that hold linear operations and the corresponding
  /// weights and biases (y = W * x + b)
  std::vector<torch::nn::Linear> _weights;
  // Number of neurons on the input layer
  const unsigned int _num_inputs;
  /// Number of neurons on the output layer
  const unsigned int _num_outputs;
  /// Hidden layer architecture
  const std::vector<unsigned int> _num_neurons_per_layer;
  /// Number of hidden layers in the neural network
  const unsigned int _num_hidden_layers;
  /// Activation functions (either one for all hidden layers or one for every layer
  /// separately)
  MultiMooseEnum _activation_function;
};

}

template <>
void dataStore<Moose::LibtorchArtificialNeuralNet>(
    std::ostream & stream,
    std::shared_ptr<Moose::LibtorchArtificialNeuralNet> & nn,
    void * context);

template <>
void dataLoad<Moose::LibtorchArtificialNeuralNet>(
    std::istream & stream,
    std::shared_ptr<Moose::LibtorchArtificialNeuralNet> & nn,
    void * context);

#endif
