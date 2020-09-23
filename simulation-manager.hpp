/*
  filename: simulation-manager.hpp

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

#ifndef STENDHAL_SIMULATION_MANAGER_HPP
#define STENDHAL_SIMULATION_MANAGER_HPP

/* Filename: simulation-manager.hpp
   Description: Implementation of the simulation-manager of stendhal package
*/


// Include C++ standard headers
#include <vector>
#include <bitset>
#include <iostream>
// Random number generators
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif // USE_PCG
#include <random>

// Include from stendhal
#include "gl_psc_exp.hpp"

namespace stendhal
{
  // Class implementation of the GL neuron model
  class SimulationManager
  {
  private:
    // initialized flag; to make sure there is only one instance
    static bool initialized;
    // Current simulation time (ms)
    static double t;
    // Simulation time step (ms)
    static double dt;
    // Vector of gl_psc_exp pointers, vector storing all created nodes
    static std::vector< class gl_psc_exp * > nodes;
    // global random number generator
#ifdef USE_PCG
    static pcg32 rng;
#else
    static std::mt19937 rng;
#endif // USE_PCG
    // global uniform random number distribution
    static std::uniform_real_distribution<> udist;
    // global poisson random number distribution
    static std::poisson_distribution<> pdist;
    // random number generator seed
    int seed{55};
  public:
    SimulationManager();
    SimulationManager(int); // Constructor with seed
    ~SimulationManager();
    // return initialized status
    bool get_initialized(void);
    // get simulation time
    double get_time(void);
    // get simulation time step
    double get_time_step(void);
    // draw uniform
    double get_uniform(void);
    // draw poisson with rate lam_
    int get_poisson(double lam_);
    // reset time
    void reset_time(void);
    // set current time to arbitrary value
    void set_time(double t_);
    // set simulation time step
    void set_dt(double delta_t);
    // set random number seed
    void set_seed(int s);
    // append node
    void append_node(class gl_psc_exp *);
    // get node from node ID
    class gl_psc_exp* get_node(int node_ID);
    // get number of nodes
    int get_n_nodes(void);
    // advance time
    void next(void);
    
  }; // class SimulationManager
} // namespace stendhal

#endif // STENDHAL_SIMULATION_MANAGER_HPP
