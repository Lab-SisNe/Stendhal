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


/* Begin Documentation

namespace: stendhal

struct stendhal::connection:
  structure containing pointer of the target neuron (gl_psc_exp class), synaptic weight and synaptic delay.
  members:
    class gl_psc_exp * target;
    double weight;
    double delay;

class stendhal::gl_psc_exp:
  class containing all the parameters, and functions of the GL neuron with
    exponentially decaying post synaptic current.

  private members:
    private variables:
      static unsigned int nodecount;
        total number of nodes created.
        nodes should be created independently for this counter to work properly
      unsigned in nodeID;
        node ID of current node
      double delta_t;
        simulation time step
        all nodes should be created with the same time step.
        there is no mechanism to check if other nodes have the same delta_t value
      std::mt19937 *prng;
        pointer for a mt19937 random number generator
      std::uniform_real_distribution<> *pudist;
        pointer for a uniform random number generator
      std::poisson_distribution<> *ppdist;
        pointer for a poisson random number generator
      bool is_local_rand;
        flag which says if the random number generators are local to the neuron
        or global to the network
      struct neuron_parameters params;
        structure containing neuron parameters.
        there is currently no function to change these parameters on the fly.
        if change in parameters is necessary, it should be done on the conde.

    state variables:
      double V_m;
        stores membrane potential for the current time step (mV)
        membrane potentials decays exponentialy to 0 with a time constant tau_m
      double I_exc;
        stores excitatory post-synaptic current for the current time step (pA)
        post-synaptic current dacays exponentialy to 0 with time constant tau_exc
      double I_inh;
        stores inhibitory post-synaptic current for the current time step (pA)
        post-synaptic current dacays exponentialy to 0 with time constant tau_inh
      double I_ext;
        stores external DC current (or stepwise) for the current time step (pA)
      int is_ref;
        refractory counter.
        is set to ref_count (see below) when the neuron fires (spikes)
        and is decreased by one for each successive time steps.
        V_m, I_exc, I_inh and I_ext are set to 0 while is_ref is positive.

    auxiliary variables:
      double rho_m;
        membrane potential leak factor: exp(-delta_t/tau_m) (unitless).
      double rho_exc;
        excitatory post-synaptic current leak factor: exp(-delta_t/tau_exc) (unitless).
      double rho_inh;
        inhibitory post-synaptic current leak factor: exp(-delta_t/tau_inh) (unitless).
      double prop_exc, prop_inh
        propagator for excitatory/inhibitory post-synaptic current (ms/pF):
        (1/C_m) * (exp(-delta_t/tau_m) - exp(-delta_t/tau_syn))/(1/tau_syn-1/tau_m)
        C_m: membrane capacitance (pF)
        tau_m: membrane potential decay time constant (ms)
        tau_syn: synaptic current decay time constant (ms). either tau_exc or tau_inh
        when this factor is multiplied by I_exc or I_syn (ms pA / pF),
        the unit conversion (A/F = V/s) results in (mV).
      double prop_step;
        propagator for stepwise external current (ms/pF):
        (1/C_m) * tau_m * (1 - exp(-delta_t/tau_m))
        similarly to prop_exc and prop_inh, when multiplied by I_ext,
        the unit conversion results in (mV)
      int ref_count = tau_ref/delta_t;
        the number of simulation steps necessary for the refractory period.

    buffer variables:
      std::vector<double> W_exc, W_inh;
        ring buffer containing excitatory (W_exc) and inhibitory (W_inh)
        synaptic weights (pA). 
        buffer position is stored in pbuff_pos and buffer size in pbuff_size
      unsigned int *pbuff_pos, *pbuff_size;
        pointer to variable containing;
          buffer position (*pbuff_pos) and
          buffer size (*pbuff_size).
        these variables can be local (independent for every neuron) or
        global (all neurons points to the same address).
        this behaviour will depend on how gl_psc_exp class is constructed, and
        stored in is_local_buff boolean variable.
        when local, pbuff_pos must be shifted independently for each neuron after
        each time step.
      bool is_local_buff;
        stores behaviour of pbuff_pos and pbuff_size.
        when true, pbuff_pos and pbuff_size are local variables.

    connection list:
      std::vector<struct connection> conn_list;
        list of connection structure (see above) containing:
	*target: pointer to post-synaptic neuron
        weight:  value of synaptic weight (pA)
        delay:  value of synaptic delay (ms)

    functions:
      double propagator_exp(tau_m, tau_syn, C_m);
        implementation of the function that calculates prop_exc and prop_inh (see above).
      double propagator_step(tau_m, C_m);
        implementation of the function that calculates prop_step (see above).
      inline double fprob(V_m);
        implementation of the firing probability function:
          V_diff = V_m - V_rheo;
          return V_diff > 0 ? pow(gamma * V_diff, r) : 0.0;

  public members:
    constructors:
      stendhal::gl_psc_exp::gl_psc_exp(seed =55);
        simplest constructor.
        node will have independent random number generators and buffer variables.
        if called without any parameter, seed is initialized to 55;
        Note that if nodes are created with the same seed, they will all behave
        with the same sequence of random numbers. That is, if they receive all the
        exact same input, they will all behave the same.
      stendhal::gl_psc_exp::gl_psc_exp(seed, *pbuff_pos, *pbuff_size);
        constructor whith global pbuff_pos and pbuff_size,
        but local random number generators.
        Note that if nodes are created with the same seed, they will all behave
        with the same sequence of random numbers. That is, if they receive all the
        exact same input, they will all behave the same.
      stendhal::gl_psc_exp::gl_psc_exp(delta_t, *prng, *pudist, *pbuff_pos, *pbuff_size);
        constructor with global random number generators and buffer variables.
        simulation time step must be the same for all neurons.
        auxiliary variables are set by calling calibrate (see below).
        State variables are reset to 0.0.
        buffer size is resized to buff_size.
        nodecount and nodeID are set.
        if class are constructed as array or vector all at once,
        all nodes in that array will have the same nodeID,
        and nodecount will not be set properly.
        * Ultimately this is the called constructor!
    destructor:
      stendhal::gl_psc_exp::~gl_psc_exp()
        deletes prng and pudist if is_local_rand
        deletes pbuff_pos and pbuff_size if is_local_buff

    public functions:
      void stendhal::gl_psc_exp::calibrate(double dt);
        update delta_t value and calculate all auxliliary variables:
          rho_m, rho_exc, rho_inh, prop_exc, prop_inh, prop_step and ref_count.
      void stendhal::gl_psc_exp::resize_buffer(unsinged int len);
        resize buffer (W_exc and W_inh) size and initialize extended elements to 0.0.
      void stendhal::gl_psc_exp::connect(stendhal::gl_psc_exp *node, double w, double d);
        append *node, w and d to the post-synaptic connection list (conn_list).
      void stendhal::gl_psc_exp::add_input(double w, unsinged int d)
        append synaptic weight (w) at position d in the ring buffer array.
          position in the ring buffer is calculated as ((*pbuff_pos)+d)%(*pbuff_size);
        if w is positive, it is added to W_exc, and if negative, to W_inh.
      void stendhal::gl_psc_exp::add_DC_input(double val);
        add val to I_ext
      double stendhal::gl_psc_exp::evaluate();
        evaluate all state variable to the next simulation step.
        return membrane potential when the neuron spiked, or 0.0 if not spiked.
        if neuron is in refractory period (is_ref > 0),
          all buffers are cleared, and no state variable are calculated.
          is_ref counter is decremented by 1.
        otherwise
          membrane potential is evolved as:
            V_m = rho_m*V_m + prop_exc*I_exc + prop_inh*I_inh + prop_step*I_ext;
          synaptic current are evolved as:
            I_exc = rho_exc*I_exc + W_exc[*pbuff_pos];
            I_inh = rho_inh*I_inh + W_inh[*pbuff_pos];
          clear input buffers: W_exc[*pbuff_pos], W_inh[*pbuff_pos], I_ext = 0.0
          Test if neuron fired:
            draw a random number (U_i) from the uniform distribution [0, 1).
            calculate firing probability (phi = fprob(V_m)).
            if (phi >= U_i)
              store V_m as the return value (ret = V_m).
              reset membrane potential and synaptic inputs:
                V_m, I_exc, I_inh = 0.0
              set refractory period counter (is_ref = ref_count).
              iterate through connection list and send spike data (function add_input)
                delay is converted from delay (ms) to timesteps as:
                  d = round(delay/delta_t)
                call (*target_node)->add_input(weight, d);
      unsingned int stendhal::gl_psc_exp::get_id();
        return the value of node_ID.
      unsigned int stendhal::gl_psc_exp::get_nodecount();
        return the value of nodecount; all nodes must have the same value.
      std::vector<struct connection>* stendhal::gl_psc_exp::get_connection();
        return pointer to the connection list vector.

 */ //End Documentation

// Include C++ standard headers
#include <vector>
#include <bitset>
#include <random> // default C++ random generator


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
    static unsigned int nodecount; // total number of nodes created;
    // node ID
    unsigned int nodeID;
    // simulation time step (ms)
    double delta_t {0.1};
    // pointer to (global) random number generator
    std::mt19937 *prng;

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
    unsigned int *pbuff_size; // pointer containing buffer size
    bool is_local_buff{false}; // flag for local or global buffer position and size
    unsigned int len; // keep a copy of buffer size

    // connection vector
    std::vector<struct connection> conn_list;
    
    // auxiliary functions
    double propagator_exp(double, double, double); // calculate propagator for exponentially decaying synaptic current
    double propagator_step(double, double); // calculate propagator for step/DC current input
    inline double frate(double); // calculate firing rate in ms
    inline double fprob(double); // calculate firing probability (unitless)
    
  public:
    // default constructor
    // step size is specified to initializa auxiliary variables. defaults to 0.1 (ms)
    // requires pointer to global random number generator engine
    // pointer to global uniform distribution, pointers to buffer position and size
    gl_psc_exp(double =0.1, std::mt19937* =NULL, std::uniform_real_distribution<>* =NULL, unsigned int* =NULL, unsigned int* =NULL);
    
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

    // Calibrate
    void calibrate(double =0.1); // calculate auxiliary variables

    // evaluate neuron dynamics in discrete time
    // uses simulation step already specified
    // if simulation time step is changed,
    // buffer data will no longer be accurate
    double evaluate(void);

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

    // Read membrane potential value
    double get_Vm();

    // Read exitatory synaptic current value
    double get_ePSC();

    // Read inhibitory synaptic current value
    double get_iPSC();

    // Read synaptic current value
    // sum of ePSC + iPSC
    double get_PSC();

    // Read extenal current value
    double get_I_ext();
    
    // Read exitatory synaptic current mediated membrane potential change
    // prop_exc * I_ext
    double get_ePSP();

    // Read inhibitory synaptic current mediated membrane potential change
    // prop_inh * I_inh
    double get_iPSP();

    // Read synaptic current mediated membrane potential changen
    // get_ePSP + get_iPSP
    double get_PSP();

    // Read external current mediated membrane potential changen
    // prop_step * I_ext
    double get_V_ext();

  }; // class gl_psc_exp

} // namespace stendhal

#endif // STENDHAL_GL_PSC_EXP_HPP
