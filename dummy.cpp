/*
  filename: dummy.cpp

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

/* Filename: dummy.cpp
   Description: Implementation of the Galves and Locherbach (1) neuron model
   in discrete time with exponentially decaying post synaptic current.

   References:
   (1) Antonio Galves and Eva Locherbach (2013).
*/


#include "dummy.hpp"

namespace stendhal
{
  unsigned int dummy::nodecount = 0;
  
  dummy::dummy()
  {
    nodecount++;
    nodeID = nodecount;
  }

  dummy::~dummy()
  {
  }

  double dummy::get_time(void)
  {
    return SimulationManager::get_time();
  }

  unsigned int dummy::get_node_ID(void)
  {
    return nodeID;
  }
  
} // namespace stendhal
