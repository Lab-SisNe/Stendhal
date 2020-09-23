/*
  filename: dummy.hpp

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

#ifndef STENDHAL_DUMMY_HPP
#define STENDHAL_DUMMY_HPP

/* Filename: dummy.hpp
   Description: Implementation of the Galves and Locherbach (1) neuron model
                in discrete time with exponentially decaying post synaptic current.

   References:
   (1) Antonio Galves and Eva Locherbach (2013).
*/


// Include C++ standard headers
#include <vector>
#include <bitset>
#include <random> // default C++ random generator

// Include PCG randon number generator if USE_PCG is defined
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif

// Include from stendhal
#include "simulation-manager.hpp"

namespace stendhal
{
  // Class implementation of the GL neuron model
  class dummy
  {
  private:
    // node count
    static unsigned int nodecount;
    // node ID
    unsigned int nodeID;
    // simulation time step (ms)
    double delta_t {0.1};
    // pointer to global random number generator
#ifdef USE_PCG
    pcg32 *prng;
#else
    std::mt19937 *prng;
#endif
    // pointer to global uniform random number distribution
    std::uniform_real_distribution<> *pudist;
    // pointer to global uniform random number distribution
    std::poisson_distribution<> *ppdist;
    // flag for local or global random number generator
    bool islocal {false};
    
    friend class SimulationManager;

  public:
    dummy(); //unsigned int*, unsigned int*); // default constructor; requires pointer to buffer position and size
    ~dummy();

    double get_time(void);
    unsigned int get_node_ID(void);
  }; // class dummy

} // namespace stendhal

#endif // STENDHAL_DUMMY_HPP
