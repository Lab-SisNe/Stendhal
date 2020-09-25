/*
  filename: gl_psc_exp.hpp

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

#ifndef STENDHAL_GL_PSC_EXP_HPP
#define STENDHAL_GL_PSC_EXP_HPP

/* Filename: gl_psc_exp.hpp
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


namespace stendhal
{
  // Define default global parameters
  const unsigned int default_buffer_length{50};

  // structure to store connection
  struct connection
  {
    class gl_psc_exp * target;
    double weight;
    double delay;
  };
    
  // Class implementation of the GL neuron model
  class gl_psc_exp
  {
  private:
    // node count
    static unsigned int nodecount;
    // node ID
    unsigned int nodeID;
    // simulation time step (ms)
    double delta_t {0.1};
    double i_delta_t {10.0};
    // pointer to (global) random number generator
#ifdef USE_PCG
    pcg32 *prng;
#else
    std::mt19937 *prng;
#endif
    // pointer to (global) uniform random number distribution
    std::uniform_real_distribution<> *pudist;
    // pointer to (global) poisson random number distribution
    std::poisson_distribution<> *ppdist;
    // flag for local or global random number generator
    bool is_local_rand {false};
    
    // define neuron params
    struct neuron_parameters
    {
      // Reset membrane potential (mV)
      double E_L {0.0};
      // Reset membrane potential (mV)
      double V_reset {0.0};
      // Membrane capacitance (pF)
      double C_m {250.0};
      // Membrane time constant (ms)
      double tau_m {10.0};
      // Time constant of excitatory postsynaptic currents (ms)
      double tau_exc {0.5};
      // Time constant of inhibitory postsynaptic currents (ms)
      double tau_inh {0.5};
      // Refractory preiod (ms)
      double tau_ref {2.0};
      // Slope of the firing probability function (1/mV)
      double gamma {0.1};
      // Curvature of the firing probability function (unitless)
      double r {0.4};
      // Rheobase potential; potential in which firing probability > 0
      double V_rheo {15.0};
    } param; // struct neuron_parameters

    // state variables
    double V_m; // membrane potential
    double I_exc; // excitatory post-synaptic current
    double I_inh; // inhibitory post-synaptic current
    double I_ext; // external input current
    int is_ref; // refractory counter

    // auxiliary variables
    double rho_m; // membrane potential leak factor
    double rho_exc; // excitatory post-synaptic current leak factor
    double rho_inh; // inhibitory post-synaptic current leak factor
    double prop_exc; // excitatory post-synaptic current propagator
    double prop_inh; // inhibitory post-synaptic current propagator
    double prop_step; // step current propagator
    int ref_count; // number of steps for refractory period
    
    // buffer
    std::vector<double> W_exc; // buffer containing excitatory synaptic weights
    std::vector<double> W_inh; // buffer containing inhibitory synaptic weights
    unsigned int *pbuff_pos; // pointer containing current position
    unsigned int *pbuff_size; // pointer containing buffe size
    bool is_local_buff{false}; // flag for local or global buffer position and size
    unsigned int len; // keep a copy of buffer size

    // connection vector
    std::vector<struct connection> conn_list;
    
    // auxiliary functions
    void calibrate(double =0.1); // calculate auxiliary variables
    double propagator_exp(double, double, double); // calculate propagator for exponentially decaying synaptic current
    double propagator_step(double, double); // calculate propagator for step/DC current input
    inline double frate(double); // calculate firing rate in ms
    /*
    inline double fprob(double, double); // calculate firing probability for arbitrary dt (in ms)
    inline double fprob(double); // calculate firing probability for dt=0.1 (ms)
    */
    
  public:
    // default constructor
    // step size is specified to initializa auxiliary variables. defaults to 0.1 (ms)
    // requires pointer to global random number generator engine
    // (pcg32 if USE_PCG is defined, mt19937 otherwise);
    // pointer to global uniform distribution, pointers to buffer position and size
#ifdef USE_PCG
    gl_psc_exp(double =0.1, pcg32* =NULL, std::uniform_real_distribution<>* =NULL, unsigned int* =NULL, unsigned int* =NULL);
#else
    gl_psc_exp(double =0.1, std::mt19937* =NULL, std::uniform_real_distribution<>* =NULL, unsigned int* =NULL, unsigned int* =NULL);
#endif
    
    // constructor with seed for local random number generators (defaults to 55)
    // note that if all nodes have the same seed, all neuron will respond the same
    // way (i.e. all random numbers generated will be the same for all neurons
    // at each time step; to avoid this, each neuron must be intialized with
    // differents seed values)
    // requires pointer to buffer position and size
    gl_psc_exp(unsigned int =55, unsigned int* =NULL, unsigned int* =NULL);

    // simplest constructor
    // all nodes will have its own variable to store
    // buffer position and size
    // when using this method, after one simulation step
    // the buffer position should be advanced in all nodes
    // similar to the above initialization method, all neurons will
    // also have its own random number generators, therefore,
    // every neuron should have different seeds to display distinct behaviour
    gl_psc_exp(unsigned int =55);
    
    // destructor
    // must delete local number generators and local buffer position if created
    ~gl_psc_exp();

    // manage buffer size
    void resize_buffer(unsigned int);

    // evaluate neuron dynamics in continuous time
    // simulation step size may be specified (defaults to 0.1)
    //double evaluate(double =0.1);

    // evaluate neuron dynamics in discrete time
    // uses simulation step already specified
    // if simulation time step is changed,
    // buffer data will no longer be accurate
    double evaluated(void);

    // return neuron ID
    unsigned int get_id(void); // return node ID

    // return node count
    unsigned int get_nodecount(void);
    
    // Create connection list;
    // append pointer to target, weight and delay
    void connect(class gl_psc_exp*, double w, double d);

    // Add synaptic input to buffer
    // delay is specified as simulation time step (discrete time mode)
    void add_input(double w, unsigned int d);

    // Add synaptic input to buffer
    // delay is specified as arrival time (continuous time mode)
    void add_input(double w, double d);

    // Add external current input
    void add_DC_input(double val);

    // Get Connection
    std::vector<struct connection>* get_connection(void);

  }; // class gl_psc_exp

} // namespace stendhal

#endif // STENDHAL_GL_PSC_EXP_HPP
