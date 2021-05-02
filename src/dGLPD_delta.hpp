/*
filename: dGLPD_delta.hpp

This code simulates the PD microcircuit model (Potjans and Diesmann, 2014)
with the GL (Galves and Locherbach, 2013) neuron in discrete time.
The leak function is an exponentially decaying function.
For the pseudo random number generation, the xoroshiro128+ algorithm was used.

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

Apr 23, 2021
*/

#ifndef STENDHAL_dGLPD_DELTA_HPP
#define STENDHAL_dGLPD_DELTA_HPP

/* Filename: dGLPD_delta.hpp
Name: PD_GL - Implementation of the Potjans and Diesmann (1) microcircuit model
              using the Galves and Locherbach (2) neuron model formalism in
              discrete time (delta_t = 0.1 ms)

References:
  (1) Tobias Potjans and Markus Diesmann (2014). 
  (2) Antonio Galves and Eva Locherbach (2013).
*/

/* Begin Documentation

namespace stendhal

class stendhal::dGLPD:
  class implementation of the Potjans and Diesmann (2014) cortical microcircuit model
  using the Galves and Locherbach (2013) stochastic neuron model formalism in discrete time.

  private members:
    variables:
      static const int N_layers {8};
        number of cortical layers
      unsigned int seed;
        seed for the global random number generator
      xoroshiro128plus prng;
        random number generator engine. uses the C++ wrapper of the
        xoroshiro128+ pseudo-random generator.
      std::uniform_real_distribution<> udist;
        real number uniform distribution generator in [0 1).
	used to draw random number to test if a neuron spikes or not.
      std::poisson_distribution<> pdist;
        poisson distribution generator
        used to draw spikes for poisson input
      std::uniform_int_distribution<> udist_int;
        integer uniform distribution generator.
        used to draw neuron indexes for creating connection.
      std::normal_distribution<> ndist;
        real normal distribution generator (mean=0.0 and sd=1.0)
        used to draw weight and delay values when creating connections.
      strut simulation_parameters sim_params;
        structure containing simulation parameters.
      struct network_parameters net_params;
        structure containing network parameters:
        unsigned int N_full[N_layers] = {20683, 5834, 21915, 5479, 4850, 1065, 14395, 2948};
          number of neurons per layer at full scale.
      double conn_prob[N_layers][N_layers]:
        connection probability matrix. columns are for source, and rows for target layers. see below for specifi values.
      unsinged int K_ext[N_layers] = {1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100};
        number of external inputs per layers.
      
        
      
 */ // End Documentation

// Include C++ standard headers
#include <vector>
#include <bitset>
#include <random> // default C++ random generator
#include <fstream>
#include <string>
#include "xoroshiro128plus.hpp" // C++ wrapper for the xoroshiro128+ PRNG

// Include from stendhal
#include "gl_psc_delta.hpp"

namespace stendhal
{
  // Class implementation of the Potjans and Diesmann (2014) cortical
  // microcircuit model using the Galves and Locherbach (2013)
  // stochastic neuron model formalism in discrete time
  class dGLPD
  {
  private:
    // Number of layers
    static const int N_layers {8};

    // simulation seed
    unsigned int seed;
    // random number generator engine
    xoroshiro128plus prng;

    // Uniform random number generator
    std::uniform_real_distribution<> udist;
    // Poisson random number generator
    std::poisson_distribution<> pdist;
    // Uniform integer distribution
    std::uniform_int_distribution<> udist_int;
    // normal distribution
    std::normal_distribution<> ndist;
    // simulation time step
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
      unsigned int K_ext[N_layers] {1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100};
      // Factor to scale indegrees
      double K_scaling {1.0};
      // Factor to scale the number of neurons
      double N_scaling {1.0};
      // Mean postsynaptic potential amplitude (in mV)
      double PSP_e {0.15};
      // Relative standard deviation of postsynaptic potential amplitude
      double PSP_sd {0.1};
      // Relative inhibitory synaptic strength
      double g {-4.0};
      // Rate of Poisonian spike generator (in Hz)
      double bg_rate {8.0};
      // Flag to turn ON or OFF poissonian background input
      // when false, use DC
      bool poisson_input {true};
      // mean delay of excitatory connections (in ms)
      double mean_delay_exc {1.5};
      // mean delay of inhibitory connections (in ms)
      double mean_delay_inh {0.75};
      // relative standard deviation of delays
      double rel_std_delay {0.5};
      // membrane capacitance
      // correct according to gl_psc_delta
      double C_m {250.0};
    } net_params; // struct network_parameters

    // buffer lenght (ms)
    unsigned int buffer_size = 50;
    unsigned int buffer_pos = 0;
    // array for scaled number of neurons per layer
    unsigned int N_scaled[N_layers];
    unsigned int N_cumsum[N_layers+1];
    // number of connection between layers
    unsigned int K_scaled[N_layers][N_layers];
    // weight matrix
    double weight_matrix[N_layers][N_layers];
    // delay matrix
    double delay_matrix[N_layers][N_layers];
    // pointer to store layer class
    //class layer_class *layer[N_layers];
    // matrix to store neuron ID (initial and final) per layer
    unsigned int pop_ID[N_layers][2];
    // vector to store pointers to gl_psc_exp class
    std::vector<gl_psc_delta*> neurons;

    // Input parameters
    double lam[N_layers];
    double DC[N_layers];

    // simulation time
    static double t;

    // output file
    std::ofstream spike_recorder;  // record spike events
    std::string spike_recorder_file {"spike_recorder.txt"};
    std::ofstream analog_recorder;  // record membrane potential and currents
    std::string analog_recorder_file {"analog_recorder.txt"};
    bool analog_rec {false};
    
  public:
    dGLPD(int =55, double =0.1, bool =false);
    ~dGLPD();

    // calculate auxiliary parameters
    void calibrate(void);
    // create neuron population
    void create_pop(void);
    // Create connection
    void connect(void);
    // Create connection from file
    void connect(std::string);
    // update buffer size
    void update_buffer_size(double);
    // Prepare method; calls calibrate, create_pop and connect(void)
    void prepare(void);
    // Simulate
    void simulate(double);

    // temporary
    double get_conn_prob(int i, int j);
    int get_N_layers(void);
    //    void check_connection(std::string);
  }; // class dGLPD
      

} // namespace stendhal

#endif //STENDHAL_dGLPD_DELTA_HPP
