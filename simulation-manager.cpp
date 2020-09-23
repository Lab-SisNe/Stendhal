/*
  filename: simulation-manager.cpp

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

/* Filename: simulation-manager.cpp
   Description: Implementation of the simulation-manager of stendhal package
*/

// Include C++ standard headers
#include <cassert>

// Include from stendhal
#include "simulation-manager.hpp"

namespace stendhal
{
  bool SimulationManager::initialized = false;
  double SimulationManager::t = 0.0;
  double SimulationManager::dt = 0.1;

  // Constructor
  SimulationManager::SimulationManager() :
    SimulationManager(55)
  {
  }
  // Constructor with seed
  SimulationManager::SimulationManager(int s) :
    seed(s)
  {
    if (initialized) {
      std::cout << "Simulation Manager already exists\nCalling destructor...";
    }
    else {
      initialized = true;
      rng.seed(seed);
    }
    delete this;
  }
  // Destructor
  SimulationManager::~SimulationManager()
  {
  }

  // return initialized status
  bool SimulationManager::get_initialized(void)
  {
    return initialized;
  }

  // return current simulation time
  double SimulationManager::get_time(void)
  {
    return t;
  }

  // return simulation time step
  double SimulationManager::get_time_step(void)
  {
    return dt;
  }

  // return uniform distribution
  double SimulationManager::get_uniform(void)
  {
    return udist(rng);
  }

  // return poisson distribution with rate lam_
  int SimulationManager::get_poisson(double lam_)
  {
    std::poisson_distribution<>::param_type p{lam_};
    return pdist(rng, p);
  }

  // reset simulation time
  void SimulationManager::reset_time(void)
  {
    t = 0.0;
  }

  // set current simulation time to arbitrary value
  void SimulationManager::set_time(double t_)
  {
    t = t_;
  }

  // set simulation time step
  void SimulationManager::set_dt(double delta_t)
  {
    dt = delta_t;
  }
  
  // set random number generator seed
  void SimulationManager::set_seed(int s)
  {
    seed = s;
    rng.seed(seed);
  }

  // append node
  void SimulationManager::append_node(class gl_psc_exp *node)
  {
    nodes.push_back(node);
  }

  // get node from node ID
  class gl_psc_exp* SimulationManager::get_node(int node_ID)
  {
    // Make shure node_ID exists
    assert ( (node_ID > 0) && (node_ID <= nodes.size()) );

    return nodes[node_ID-1];
  }

  // get number of nodes
  int SimulationManager::get_n_nodes(void)
  {
    return nodes.size();
  }

  // advance time
  void SimulationManager::next(void)
  {
    t += dt;
  }
  
} // namespace stendhal

