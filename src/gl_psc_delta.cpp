/*
  filename: gl_psc_delta.cpp

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

/* Filename: gl_psc_delta.cpp
   Description: Implementation of the Galves and Locherbach (1) neuron model
   in discrete time with delta shape post synaptic current.

   References:
   (1) Antonio Galves and Eva Locherbach (2013).
*/

#include <iostream>
#include "gl_psc_delta.hpp"

namespace stendhal
{
  unsigned int gl_psc_delta::nodecount = 0;
  
  // Constructor with delta_t, len and global random number genertors
  // default constructor
  // step size is specified to initializa auxiliary variables. defaults to 0.1 (ms)
  // requires pointer to global random number generator engine
  // pointer to global uniform distribution, pointers to buffer position and size
  gl_psc_delta::gl_psc_delta(double dt, xoroshiro128plus *rng, std::uniform_real_distribution<> *dist, unsigned int *pos, unsigned int *size)
  :
    pbuff_pos(pos), // pointer to buffer position variable
    pbuff_size(size), // pointer to buffer size variable
    prng(rng), // pointer to global randon number generator engine
    pudist(dist) // pointer to uniform distribution generator
  {
    // calculate auxiliary variables that depends on dt
    calibrate(dt);
    // pointer to global random number generator engine
    prng = rng;
    // pointer to random number distribution
    pudist = dist;
    // Reset state variables
    V_m = 0.0;
    I_exc = 0.0;
    I_inh = 0.0;
    I_ext = 0.0;
    is_ref = 0;
    // set buffer size
    resize_buffer(*size);
    // increment total number of nodes created (static variable)
    nodecount++;
    // record current node ID
    nodeID = nodecount;
  }

  // constructor with seed for local random number generators (defaults to 55)
  // note that if all nodes have the same seed, all neuron will respond the same
  // way (i.e. all random numbers generated will be the same for all neurons
  // at each time step; to avoid this, each neuron must be intialized with
  // differents seed values)
  // requires pointer to buffer position and size
  gl_psc_delta::gl_psc_delta(unsigned int seed, unsigned int *pos, unsigned int *size) :
    gl_psc_delta(0.1, new xoroshiro128plus, new std::uniform_real_distribution<>, pos, size)
  {
    // set flag for local random numbers
    is_local_rand = true;
  }

  // simplest constructor
  // all nodes will have its own variable to store
  // buffer position and size
  // when using this method, after one simulation step
  // the buffer position should be advanced in all nodes
  // similar to the above initialization method, all neurons will
  // also have its own random number generators, therefore,
  // every neuron should have different seeds to display distinct behaviour
  gl_psc_delta::gl_psc_delta(unsigned int seed) :
    gl_psc_delta(seed, new unsigned int, new unsigned int)
  {
    (*pbuff_pos) = 0.0;
    (*pbuff_size) = default_buffer_length;
    is_local_buff = true;
  }

  // destructor
  gl_psc_delta::~gl_psc_delta()
  {
    // delete random number objects if local
    if (is_local_rand) {
      delete prng;
      delete pudist;
    }
    // delete buffer position and size objects if local
    if (is_local_buff) {
      delete pbuff_pos;
      delete pbuff_size;
    }
    // nothing need to be done if objects are global
  }

  // calibrate by specifing dt
  void gl_psc_delta::calibrate(double dt)
  {
    delta_t = dt; // update dt value

    // calculate membrane potential leak factor
    rho_m = std::exp(-delta_t/param.tau_m);
    // calculate DC propagator
    prop_step = propagator_step(param.tau_m, param.C_m);
    // number of steps for refractory period
    ref_count = (int)(param.tau_ref/delta_t);
  } // calibrate

  // step/gDC propagator
  double gl_psc_delta::propagator_step(double tau_m, double C_m)
  {
    return (1/C_m)*tau_m*(1 - std::exp(-delta_t/tau_m));
  } // step propagator

  // Firing rate in spikes/ms
  inline double gl_psc_delta::frate(double V_m)
  {
    double V_diff = V_m - param.V_rheo;
    if (V_diff > 0)
      return 10.0*std::pow(param.gamma * V_diff, param.r);
    else
      return 0.0;
  } // phi

  // Firing probability (unitless)
  inline double gl_psc_delta::fprob(double V_m)
  {
    double V_diff = V_m - param.V_rheo;
    if (V_diff > 0)
      return std::pow(param.gamma * V_diff, param.r);
    else
      return 0.0;
  } // phi

  // resize buffer
  void gl_psc_delta::resize_buffer(unsigned int len)
  {
    W_exc.resize(len, 0.0);
    W_inh.resize(len, 0.0);
  } // resize_buffer

  // evaluate with pre-defined step size (discrete algorithm)
  double gl_psc_delta::evaluate(void)
  {
    double ret {0.0};
    if (is_ref>0) {
      is_ref--;
      // clear input buffer
      W_exc[*pbuff_pos] = 0.0;
      W_inh[*pbuff_pos] = 0.0;
      I_ext = 0.0;
    }
    else {
      // Evolve synaptic currents;
      I_exc = W_exc[*pbuff_pos];
      I_inh = W_inh[*pbuff_pos];
      // Evolve membrane potential to the next simulation step
      V_m = rho_m*V_m + I_exc + I_inh + prop_step*I_ext;
      // Clear input buffer
      W_exc[*pbuff_pos] = 0.0;
      W_inh[*pbuff_pos] = 0.0;
      //I_ext = 0.0;
      // Draw random number to test if neuron fires or not
      // note that udist and rng are pointers
      double U_i = (*pudist)(*prng);
      // Calculate firing rate for the time window in question
      // This could be taken as the rate of a poisson process
      // Could we just draw a poisson with rate phi to obtain
      // the number of spikes for the time window delta_t?
      // double phi = delta_t*frate(V_m);
      double phi = fprob(V_m);
      // Neuron fired!
      if (phi >= U_i) {
	// Store membrane potential when neuron fired
	ret = V_m;
	// Reset membrane potential
	V_m = param.V_reset;
	// clear synaptic current
	I_exc = 0.0;
	I_inh = 0.0;
	// set refractory period counter
	is_ref = ref_count;
	// Send spike
	// iterate through all connection list
	for (std::vector<struct connection>::iterator it=conn_list.begin(); it!=conn_list.end(); it++) {
	  // call add_input method of the post-synaptic neuron
	  (*it).target->add_input((*it).weight, (*it).d);
	}
	// return potential in which the neuron fired
	return ret;
      } // neuron fired
    } // not refractory
    // return zero when neuron did not fire
    return ret;
  } // evaluate pre-defined step size (discrete algorithm)

  /*
  // evaluate with specific step size (continuous algorithm)
  double gl_psc_delta::evaluate(double dt)
  {
    double ret {0.0};
    
    calibrate(dt); // update auxiliary variables to the new delta_t

    // using same algorithm as discrete time for the moment. need update
    return evaluate();
  } // evaluate
  */
  
  // return node ID
  unsigned int gl_psc_delta::get_id(void)
  {
    return nodeID;
  }

  // return node count
  unsigned int gl_psc_delta::get_nodecount(void)
  {
    return nodecount;
  }
    
  // Create connection list
  void gl_psc_delta::connect(class gl_psc_delta* node, double w, double d, int pos)
  {
    // append post neuron pointer, weight (pA) and delay (ms)
    conn_list.push_back({node, w, d, pos});
  }

  // Receive input (discrete time)
  void gl_psc_delta::add_input(double w, unsigned int d)
  {
    //assert (d<(*pbuff_size));
    // calculate position on the ring buffer
    // modulo of (current_pos + delay) / buffer_length
    unsigned int pos = ((*pbuff_pos)+d)%(*pbuff_size);
    if (w>0) {
      W_exc[pos] += w; // add weight to W_exc if w is positive
    }
    else {
      W_inh[pos] += w; // add weight to W_inh if w is negative
    }
  }

  // add External current input
  void gl_psc_delta::add_DC_input(double val)
  {
    I_ext = val;
  }
  
  // Receive input (continuous time)
  void gl_psc_delta::add_input(double w, double d)
  {
    // ToDo: Implement code for continuous time
    // This code may not be exactly here
  }

  std::vector<struct connection>* gl_psc_delta::get_connection(void) {
    return &conn_list;
  }

  // Read membrane potential value
  double gl_psc_delta::get_Vm()
  {
    return V_m;
  }

  // Read exitatory synaptic current value
  double gl_psc_delta::get_ePSC()
  {
    return I_exc;
  }

  // Read inhibitory synaptic current value
  double gl_psc_delta::get_iPSC()
  {
    return I_inh;
  }

  // Read synaptic current value
  // sum of I_exc + I_inh
  double gl_psc_delta::get_PSC()
  {
    return I_exc + I_inh;
  }

  // Read extenal current value
  double gl_psc_delta::get_I_ext()
  {
    return I_ext;
  }

  // Read exitatory synaptic current mediated membrane potential change
  // prop_exc * I_exc
  double gl_psc_delta::get_ePSP()
  {
    return I_exc;
  }

  // Read exitatory synaptic current mediated membrane potential change
  // prop_inh * I_inh
  double gl_psc_delta::get_iPSP()
  {
    return I_inh;
  }

  // Read exitatory synaptic current mediated membrane potential changen
  // prop V_exc + V_inh
  double gl_psc_delta::get_PSP()
  {
    return I_exc + I_inh;
  }

  // Read external current mediated membrane potential change
  // prop_step * I_ext
  double gl_psc_delta::get_V_ext()
  {
    return prop_step * I_ext;
  }

  // Get membrane capacitance
  double gl_psc_delta::get_C_m()
  {
    return param.C_m;
  }

} // namespace stendhal
