/*
filename: network.hpp

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

#ifndef STENDHAL_NETWORK_HPP
#define STENDHAL_NETWORK_HPP

/* Filename: network.hpp
Name: network - Implementation of neuron network using the
                Galves and Locherbach (1) neuron model.

References:
  (1) Antonio Galves and Eva Locherbach (2013).
*/

#include <vector>

namespace stendhal
{
  class network
  {
  public:
    // Store number of nodes
    static int n
    network() { n++; }

    std::vector target;
    std::vector weight;
    std::vector delay;
  };

  class node
  {
  public:
    // Store number of nodes
    static int n
    node() { n++; }
  }
} // namespace stendhal

#endif STENDHAL_NETWORK_HPP
