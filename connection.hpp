/*
  filename: connection.hpp

  This code simulates the PD microcircuit model (Potjans and Diesmann, 2014)
  with the GL (Galves and Locherbach, 2013) neuron in discrete time.

  This code is part of STENDHAL package.
  A NeuroMat package for simulate neuron networks, and analyse data.

  Director:
  Antonio Galves

  Developers:
  Nilton L. Kamiji
  Christophe Pouzat
  Renan Shimoura

  Contributors:
  Karine Guimaraes
  Aline Duarte
  Jorge Stolfi
  Antonio Roque

  Aug 11, 2020
*/

#ifndef STENDHAL_CONNECTION_HPP
#define STENDHAL_CONNECTION_HPP

/* Filename: connection.hpp
   Description - Implementation of synaptic connection between nodes.
*/


// Include C++ standard headers
#include <vector>
#include <bitset>
#include <random> // default C++ random generator

namespace stendhal
{
  // Class implementation synapse
  class synapse
  {
  private:
    // presynaptic node ID
    unsigned int pre;
    // postsynaptic node ID
    unsigned int post;
    // synaptic weight
    double weight;
    // synaptic delay
    double delay;
    
    // buffer
    std::vector<bool> X; // buffer containing spike train
    unsigned int *buff_pos; // pointer containing current position
    unsigned int *buff_size; // pointer containing buffer size

  public:
    synapse(); // default constructor
    ~synapse();
    void resize_buffer(int); // resize buffer length
  }; // class synapse

} // namespace stendhal

#endif // STENDHAL_CONNECTION_HPP
