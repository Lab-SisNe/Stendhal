/*
  filename: gl_psc_exp.cpp

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

/* Filename: gl_psc_exp.cpp
   Description: Implementation of the Galves and Locherbach (1) neuron model
   in discrete time with exponentially decaying post synaptic current.

   References:
   (1) Antonio Galves and Eva Locherbach (2013).
*/


#include "gl_psc_exp.hpp"

namespace stendhal
{
  unsigned int gl_psc_exp::nodecount = 0;
  
  // Constructor with delta_t, len and global random number genertors
  // default constructor
  // step size is specified to initializa auxiliary variables. defaults to 0.1 (ms)
  // requires pointer to global random number generator engine
  // (pcg32 if USE_PCG is defined, mt19937 otherwise);
  // pointer to global uniform distribution, pointers to buffer position and size
#ifdef USE_PCG
  gl_psc_exp::gl_psc_exp(double dt, pcg32 *rng, std::uniform_real_distribution<> *dist, unsigned int *pos, unsigned int *size)
#else
  gl_psc_exp::gl_psc_exp(double dt, std::mt19937 *rng, std::uniform_real_distribution<> *dist, unsigned int *pos, unsigned int *size)
#endif
  :
    pbuff_pos(pos),
    pbuff_size(size)
  {
    if (dt != delta_t) {
      delta_t = dt;
      i_delta_t = 1.0/dt; // multiplication is faster than division; used when divion by dt is necessary
    }
    if (len != *size) {
      len = *size; // keep a copy of global buffer length; may not be necessary
    }
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
  gl_psc_exp::gl_psc_exp(unsigned int seed, unsigned int *pos, unsigned int *size) :
    pbuff_pos(pos),
    pbuff_size(size)
  {
    // create local random number generator engine
#ifdef USE_PCG
    prng = new pcg32(seed);
#else
    prng = new std::mt19937(seed);
#endif
    // create local random number distribution
    pudist = new std::uniform_real_distribution<>;
    // set flag for local random numbers
    is_local_rand = true;
    // set buffer size
    resize_buffer(*size);
    // increment total number of nodes created (static variable)
    nodecount++;
    // record current node ID
    nodeID = nodecount;
  }

  // simplest constructor
  // all nodes will have its own variable to store
  // buffer position and size
  // when using this method, after one simulation step
  // the buffer position should be advanced in all nodes
  // similar to the above initialization method, all neurons will
  // also have its own random number generators, therefore,
  // every neuron should have different seeds to display distinct behaviour
  gl_psc_exp::gl_psc_exp(unsigned int seed) :
    gl_psc_exp(seed, new unsigned int, new unsigned int)
  {
    (*pbuff_pos) = 0.0;
    (*pbuff_size) = default_buffer_length;
    is_local_buff = true;
  }

  // destructor
  gl_psc_exp::~gl_psc_exp()
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
  void gl_psc_exp::calibrate(double dt)
  {
    delta_t = dt; // update dt value
    i_delta_t = 1.0/dt;

    // calculate membrane potential leak factor
    rho_m = std::exp(-delta_t/param.C_m);
    // calculate excitatory post-synaptic current leak factor
    rho_exc = std::exp(-delta_t/param.tau_exc);
    // calculate inhibitory post-synaptic current leak factor
    rho_inh = std::exp(-delta_t/param.tau_inh);
    // calculate excitatory post-synaptic current propagator
    prop_exc = propagator_exp(param.tau_m, param.tau_exc, param.C_m);
    // calculate inhibitory post-synaptic current propagator
    prop_inh = propagator_exp(param.tau_m, param.tau_inh, param.C_m);
    // calculate DC propagator
    prop_step = propagator_step(param.tau_m, param.C_m);
    // number of steps for refractory period
    ref_count = (int)(param.tau_ref*i_delta_t);
  } // calibrate

  // exponential propagator
  double gl_psc_exp::propagator_exp(double tau_m, double tau_syn, double C_m)
  {
    if (tau_m == tau_syn)
      return delta_t * std::exp(-delta_t/tau_m) / C_m;
    else
      return (1/C_m)*(tau_m*tau_syn/(tau_m-tau_syn))*(std::exp(-delta_t/tau_m) - std::exp(-delta_t/tau_syn));
  } // exponential propagator

  // step/DC propagator
  double gl_psc_exp::propagator_step(double tau_m, double C_m)
  {
    return (1/C_m)*tau_m*(1 - std::exp(-delta_t/tau_m));
  } // step propagator

  // Firing rate in spikes/ms
  inline double gl_psc_exp::frate(double V_m)
  {
    double V_diff = V_m - param.V_rheo;
    if (V_diff > 0)
      return 10.0*std::pow(param.gamma * V_diff, param.r);
    else
      return 0.0;
  } // phi

  /*
  // Firing probability for arbitrary dt (in ms)
  inline double gl_psc_exp::fprob(double V_m, double d_t)
  {
    return d_t*frate(V_m);
  } // fprob
  
  // Firing probability for dt=0.1 ms
  inline double gl_psc_exp::fprob(double V_m)
  {
    return 0.1*frate(V_m);
  } // phi
  */
  
  // resize buffer
  void gl_psc_exp::resize_buffer(unsigned int len)
  {
    W_exc.resize(len, 0.0);
    W_inh.resize(len, 0.0);
  } // resize_buffer

  // evaluate with pre-defined step size (discrete algorithm)
  double gl_psc_exp::evaluated(void)
  {
    double ret {0.0};
    if (is_ref>0)
      is_ref--;
    else {
      // Evolve membrane potential to the next simulation step
      V_m = rho_m*V_m + prop_exc*I_exc + prop_inh*I_inh + prop_step*I_ext;
      // Evolve synaptic currents;
      I_exc = rho_exc*I_exc + W_exc[*pbuff_pos];
      I_inh = rho_inh*I_inh + W_inh[*pbuff_pos];
      // Clear input buffer
      W_exc[*pbuff_pos] = 0.0;
      W_inh[*pbuff_pos] = 0.0;
      I_ext = 0.0;
      // Draw random number to test if neuron fires or not
      // note that udist and rng are pointers
      double U_i = (*pudist)(*prng);
      // Calculate firing rate for the time window in question
      // This could be taken as the rate of a poisson process
      // Could we just draw a poisson with rate phi to obtain
      // the number of spikes for the time window delta_t?
      double phi = delta_t*frate(V_m);
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
	for (std::vector<connection>::iterator it=conn_list.begin(); it!=conn_list.end(); it++) {
	  // Convert delay to number of steps
	  // increment one as neuron spiked at time t+dt
	  // and buffer position is still at time t
	  // do not round, as spikes arriving with a delay of 0.19
	  // should have the same effect as a delay of 0.1,
	  // and not as 0.2 or 0.21
	  //unsigned int d = std::round((*it).delay*i_delta_t)+1;
	  unsigned int d = (int)((*it).delay*i_delta_t)+1;
	  // call add_input method of the post-synaptic neuron
	  (*it).target->add_input((*it).weight, d);
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
  double gl_psc_exp::evaluate(double dt)
  {
    double ret {0.0};
    
    calibrate(dt); // update auxiliary variables to the new delta_t

    // using same algorithm as discrete time for the moment. need update
    return evaluate();
  } // evaluate
  */
  
  // return node ID
  unsigned int gl_psc_exp::get_id(void)
  {
    return nodeID;
  }

  // return node count
  unsigned int gl_psc_exp::get_nodecount(void)
  {
    return nodecount;
  }
    
  // Create connection list
  void gl_psc_exp::connect(class gl_psc_exp* t, double w, double d)
  {
    // append post neuron pointer, weight (pA) and delay (ms)
    conn_list.push_back({t, w, d});
  }

  // Receive input (discrete time)
  void gl_psc_exp::add_input(double w, unsigned int d)
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
  void gl_psc_exp::add_DC_input(double val)
  {
    I_ext += val;
  }
  
  // Receive input (continuous time)
  void gl_psc_exp::add_input(double w, double d)
  {
    // ToDo: Implement code for continuous time
    // This code may not be exactly here
  }

} // namespace stendhal
