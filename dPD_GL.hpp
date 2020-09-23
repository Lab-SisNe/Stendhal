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
#include <fstream>

// Include pcg if USE_PCG is defined; compile with -DUSE_PCG
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp" // PCG random number generator
#endif

// Include from stendhal
#include "gl_psc_exp.hpp"

namespace stendhal
{
  // Class implementation of the Potjans and Diesmann (2014) cortical
  // microcircuit model using the Galves and Locherbach (2013)
  // stochastic neuron model formalism in discrete time
  class dPD_GL
  {
  private:
    // Number of layers
    static const int N_layers {8};

    // simulation seed
    unsigned int seed;
    // random number generator engine
#ifdef USE_PCG
    pcg32 rng;
#else
    std::mt19937 rng;
#endif
    // Uniform random number generator
    std::uniform_real_distribution<> udist;
    // Poisson random number generator
    std::poisson_distribution<> pdist;
    // Uniform integer distribution
    std::uniform_int_distribution<> udist_int;
    // normal distribution
    std::normal_distribution<> ndist;
    // simulatio time step
    double delta_t;
    
    // define simulation parameters
    struct simulation_parameters
    {
      // simulation time step (ms)
      double delta_t {0.1};
      // simulation time (ms)
      double t_sim {1000.0};
      // random number generator seed
      unsigned int seed {55};
      
      // number of threds
      // int n_threads {1};
      // number of MPI processes; 
      // int n_mpi {1};
    } sim_params; // struct simulation_parameters

    // define network params
    struct network_parameters
    {
      // Number of neurons per layer
      // L23E, L23I, L4E, L4I, L5E, L5I, L6E, L6I
      unsigned int N_full[N_layers] {20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948};
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
      unsigned int K_ext[N_layers] {1600, 1500, 2100, 1900, 2000, 1200, 2900, 2100};
      // Factor to scale indegrees
      double K_scaling {1.0};
      // Factor to scale the number of neurons
      double N_scaling {1.0};
      // Mean postsynaptic current amplitude (in pA)
      double PSP_e {87.8};
      // Relative standard deviation of postsynaptic current amplitude
      double PSP_sd {0.1};
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
      // synaptic decay time constant
      // Correct accordingly to parameter used in gl_psc_exp
      double tau_syn_ex {0.5};
    } net_params; // struct network_parameters

    // buffer lenght (ms)
    double buffer_len_ms = 5.0;
    unsigned int buffer_size = 50;
    unsigned int buffer_pos = 0;
    // array for scaled number of neurons per layer
    unsigned int N_scaled[N_layers];
    unsigned int N_cumsum[N_layers+1];
    // number of connection between layers
    unsigned int K_scaled[N_layers][N_layers];
    // pointer to store layer class
    class layer_class *layer[N_layers];
    // matrix to store neuron ID (initial and final) per layer
    unsigned int pop_ID[N_layers][2];
    // vector to store pointers to gl_psc_exp class
    std::vector<gl_psc_exp*> neurons;

    // Input parameters
    double lam[N_layers];
    double DC[N_layers];

    // simulation time
    static double t;

    // output file
    std::ofstream spike_recorder;
    std::string spike_recorder_file {"spike_recorder.txt"};
    
  public:
    dPD_GL(int =55);
    ~dPD_GL();

    // calculate auxiliary parameters
    void calibrate(void);
    // create neuron population
    void create_pop(void);
    // Create connection
    void connect(void);
    // Simulate
    void simulate(double);
  }; // class dPD_GL
      

} // namespace stendhal

#endif //STENDHAL_dPD_GL_HPP
