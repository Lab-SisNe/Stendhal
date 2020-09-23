/*
filename: PD_GL.hpp

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

#ifndef STENDHAL_dPD_GL_HPP
#define STENDHAL_dPD_GL_HPP

/* Filename: dPD_GL.hpp
Name: PD_GL - Implementation of the Potjans and Diesmann (1) microcircuit model
              using the Galves and Locherbach (2) neuron model formalism in
              discrete time (delta_t = 0.1 ms)

References:
  (1) Tobias Potjans and Markus Diesmann (2014). 
  (2) Antonio Galves and Eva Locherbach (2013).
*/

// Flag to select random number generator
// comment out to use PCG
//#define USE_PCG

// Include C++ standard headers
#include <vector>
#include <bitset>
#include <random> // default C++ random generator
#ifdef USE_PCG
#include "pcg_cpp/pcg_random.hpp" // PCG random number generator
#endif

namespace stendhal
{
  // Class implementation of the Potjans and Diesmann (2014) cortical
  // microcircuit model using the Galves and Locherbach (2013)
  // stochastic neuron model formalism in discrete time
  class dPD_GL
  {
  private:
    // define simulation parameters
    struct simulation_parameters
    {
      // simulation time step (ms)
      double delta_t {0.1};
      // simulation time (ms)
      double t_sim {1000.0};
      // random number generator seed
      int seed {55};
      // random number generator
#ifdef USE_PCG
      pcg32 rng;
#else
      std::mt19937 rng;
#endif
      // number of threds
      // int n_threads {1};
      // number of MPI processes; 
      // int n_mpi {1};
    }; // struct simulation_parameters

    // define network params
    struct network_parameters
    {
      // Number of layers
      int N_layers {8};
      // Number of neurons per layer
      // L23E, L23I, L4E, L4I, L5E, L5I, L6E, L6I
      int N_full[8] {20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948};
      // Connection probability
      //                           from
      //            L23E L23I L4E L4I L5E L5I L6E L6I
      //     L23E
      //     L23I
      //     L4E
      //     L4I
      // to  L5E
      //     L5I
      //     L6E
      //     L6I
      double conn_prob[N_layers][N_layers] {
	{0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0., 0.0076, 0.},
	{0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0., 0.0042, 0.},
	{0.0077, 0.0059, 0.0497, 0.135, 0.0067, 0.0003, 0.0453, 0.},
	{0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0., 0.1057, 0.},
	{0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.},
	{0.0548, 0.0269, 0.0257, 0.0022, 0.06, 0.3158, 0.0086, 0.},
	{0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252},
	{0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}
	};
      // Number of external connections to the different populations
      double K_ext[N_layers] {1600, 1500, 2100, 1900, 2000, 1200, 2900, 2100};
      // Factor to scale indegrees
      double K_scaling {1.0};
      // Factor to scale the number of neurons
      double N_scaling {0.1};
      // Mean postsynaptic current amplitude (in pA)
      double PSP_e {87.8};
      // Relative inhibitory synaptic strength
      double g {-4.0};
      // Rate of Poisonian spike generator (in Hz)
      double bg_rate {8.0};
      // Flag to turn ON or OFF poissonian background input
      // when false, use DC
      bool poisson_input {false};
      // mean delay of excitatory connections (in ms)
      double mean_delay_exc {1.5};
      // mean delay of inhibitory connections (in ms)
      double mean_delay_inh {0.8};
      // relative standard deviation of delays
      double rel_std_delay {0.5};
    }; // struct network_parameters

    // buffer lenght (ms)
    double buffer_len_ms = 5.0;
    // array for scaled number of neurons per layer
    int N_scaled[N_layers];
    // pointer to store layer class
    class layer_class *layer[N_layers];

    // Define layer class
    class layer_class
    {
    private:
      // define neuron params
      struct neuron_parameters
      {
	// mean initial membrane potential (mV)
	double V0_mean {7.0};
	// standard deviation of the initial membrane potential (V)
	double V0_sd {10.0};
	// Reset membrane potential (mV)
	double E_L {0.0};
	// LIF Threshold potential (mV)
	double V_th {15.0};
	// Reset membrane potential (mV)
	double V_reset {0.0};
	// Membrane capacitance (pF)
	double C_m {250.0};
	// Membrane time constant (ms)
	double tau_m {10.0};
	// Time constant of excitatory postsynaptic currents (ms)
	double tau_syn_ex {0.5};
	// Time constant of inhibitory postsynaptic currents (ms)
	double tau_syn_in {0.5};
	// Refractory preiod (ms)
	double tau_ref {2.0};
	// Slope of the firing probability function (1/mV)
	double gamma {0.1};
	// Curvature of the firing probability function (unitless)
	double r {0.4};
	// Rheobase potential; potential in which firing probability > 0
	double V_rheo {15.0};
      }; // struct neuron_parameters

      // number of elements
      int size;
      // state variables
      double *V_m; // membrane potential
      double *I_syn_ex; // excitatory post-synaptic current
      double *I_syn_in; // inhibitory post-synaptic current
      double *I_ext; // external input current
      int *is_ref; // refractory counter

      // buffer
      std::vector<double> *W_ex;
      std::vector<double> *W_in;
      
    public:
      layer_class(int N, int len);
      /*
	size(N)
      {
	V_m = new double[size];
	I_syn_ex = new double[size];
	I_syn_in = new double[size];
	I_ext = new double[size];
	is_ref = new int[size];
	W_ex = new std::vector<double>[size];
	W_in = new std::vector<double>[size];
	for (int i=0; i<size; i++) {
	  W_ex[i].resize(len, 0.0);
	  W_in[i].resize(len, 0.0);
	}
      };
      */
      ~layer_class();
      /*
      {
	delete []V_m;
	delete []I_syn_ex;
	delete []I_syn_in;
	delete []I_ext;
	delete []is_ref;
	delete []W_ex;
	delete []W_in;
      };
      */

      void resize_buffer(int len)
      {
	for (int i=0; i<size; i++) {
	  W_ex[i].resize(len, 0.0);
	  W_in[i].resize(len, 0.0);
	}
      };
    }; // class layer_class

    // Define synapse class
    class synapse_class
    {
    private:
      std::bitset<128> input_buffer; // Create a 128 bit buffer to store spike event
      unsigned int delay; // synaptic transmission delay in time steps
      double weight; // synaptic weith
      
    public:
      synapse_class();
      ~synapse_class();
    }; // class synapese_class
      
    // Define external input class
    class external_input
    {
    private:
      //int size; // number of elements
      bool is_poisson; // flag 
      int poisson_stim; // mean spike count per time step for poisson input
      double DC_stim; //
      //int *target; // array containing target neuron ID.
      std::function<void(void)> apply_stim;  
    public:
      external_input(int N);
      ~external_input();
      void poisson_input(void);
      void DC_input(void);
    }; // class external_input
    
  public:
    dPD_GL();
    dPD_GL(int);
    ~dPD_GL();

    void print_N_full(void);

    void print_conn_prob(int, int);
  }; // class dPD_GL

      

} // namespace stendhal

#endif //STENDHAL_dPD_GL_HPP
