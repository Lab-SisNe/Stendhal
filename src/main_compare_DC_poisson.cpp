/*
filename: main_compare_DC_poisson.cpp

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

// Include C++ standard headers
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <chrono>

// Include stendhal headers
#include "stendhalconfig.h"
#include "gl_psc_exp.hpp"

/* Filename: main_dPD_GL.cpp
main function

*/

// Usage message
static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name << " seed t_sim <conn_csv_data>" << std::endl
	    << "  seed: seed for random number generator (integer)" << std::endl
	    << "  t_sim: simulation time in (ms)" << std::endl
	    << "  delta_t: simulation time step (ms)" << std::endl
	    << "  PSC_e: synaptic weight (pA)" << std::endl    
    //	    << "  DC_start: start value of DC current:" << std::endl
    //	    << "  DC_end: final value of DC current:" << std::endl
    //	    << "  DC_step: step value of DC current:" << std::endl
	    << std::endl
	    << " Simulate N neurons to constant DC current input."
	    << std::endl;
}

int main ( int argc, char* argv[] )
{
  // Get input arguments
  if ((argc != 5)) {
    show_usage(argv[0]);
    return 1;
  }
  int seed = std::atoi(argv[1]);
  double t_sim = std::atof(argv[2]);
  double delta_t = std::atof(argv[3]);
  double PSC_e = std::atof(argv[4]);
  //  double DC_start = std::atof(argv[4]);
  //  double DC_end = std::atof(argv[5]);
  //  double DC_step = std::atof(argv[6]);
  std::cout << "Random Number Generator Engine: ";
  std::cout << "mt19937" << std::endl;
  std::cout << "seed: " << seed
	    << ", Simulation time: "
	    << t_sim << " ms"
	    << ", Simulation time step: "
	    << delta_t << " ms" << std::endl;
  //  std::cout << "DC_start = " << DC_start
  //	    << ", DC_end = " << DC_end
  //	    << ", DC_step = " << DC_step << std::endl;

  // neurons, Kext, N
#define N_LAYER 8
  std::vector<stendhal::gl_psc_exp*> neurons;
  double DC[N_LAYER], lam[N_LAYER];
  unsigned int K_ext[8] {1600, 1500, 2100, 1900, 2000, 1900, 2900, 2100};
  unsigned int buffer_pos {0};
  unsigned int buffer_size {5};
  double t = 0.0;
  double bg_rate {8.0}; // Hz
  double PSC_e {87.8}; // pA
  double tau_syn {0.5}; // ms
  
  // random number generators
  std::mt19937 rng;
  // uniform random number generator
  std::uniform_real_distribution<> udist;
  // Poisson random number generator
  std::poisson_distribution<> pdist;
  // parameter variables for poisson distribution
  std::poisson_distribution<>::param_type pparam;
  // number of poissonian inputs
  unsigned int p;

  rng.seed(seed);
  
  // Timer
  //  std::chrono::time_point<std::chrono::steady_clock> start;
  //  std::chrono::time_point<std::chrono::steady_clock> end;
  //  std::chrono::duration<double> diff;

  // outpuf file
  std::ofstream DC_recorder;
  std::string DC_recorder_file {"DC_recorder.txt"};
  std::ofstream spike_recorder;
  std::string spike_recorder_file {"spike_recorder.txt"};

  // open output file
  DC_recorder.open(DC_recorder_file);
  DC_recorder << "# dt = " << delta_t << std::endl;
  spike_recorder.open(spike_recorder_file);
  spike_recorder << "# dt = " << delta_t << std::endl;
  spike_recorder << "# time,  neuron_ID,  V_m" << std::endl;
  
  // Create neurons for DC input
  for (int i=0; i<N_LAYER; i++) {
    // append new gl_psc_exp instance
    neurons.push_back(new stendhal::gl_psc_exp(delta_t, &rng, &udist, &buffer_pos, &buffer_size));
    // DC amplitude per layer
    DC[i] = bg_rate * K_ext[i] * PSC_e * tau_syn * 1e-3;
    DC_recorder << DC[i] << std::endl;
  }
  DC_recorder << "# lambda for poisson input:" << std::endl;
  // Create neurons for poisson input
  for (int i=0; i<N_LAYER; i++) {
    // append new gl_psc_exp instance
    neurons.push_back(new stendhal::gl_psc_exp(delta_t, &rng, &udist, &buffer_pos, &buffer_size));
    // lambda for poisson distribution per layer
    lam[i] = bg_rate * K_ext[i] * delta_t * 1e-3;
    DC_recorder << lam[i] << std::endl;
  }

  // simulate
  double V_spiked;
  while (t<=t_sim) {
    t += delta_t;
    for (int i=0; i<N_LAYER; i++) {
      // Simulate DC neurons
      neurons[i]->add_DC_input(DC[i]); // add DC input in neuron i
      V_spiked = neurons[i]->evaluated(); // evaluate neuron i
      // store data if neuron[i] fired
      if (V_spiked > 0) {
	spike_recorder << (int)std::round(t/delta_t) << ", "
		       << neurons[i]->get_id() << ", "
		       << V_spiked << std::endl;
      }
      // Simulate poisson neuron
      pparam = std::poisson_distribution<>::param_type (lam[i]); // set param
      p = pdist(rng, pparam); // draw poisson distribution with rate lam[i]
      neurons[N_LAYER + i]->add_input(p*PSC_e,(unsigned int)0); // set poisson input in enruon N_LAYER + i
      V_spiked = neurons[N_LAYER + i]->evaluated(); // evaluate neuron N_LAYER + i
      if (V_spiked > 0) {
	spike_recorder << (int)std::round(t/delta_t) << ", "
		       << neurons[N_LAYER + i]->get_id() << ", "
		       << V_spiked << std::endl;
      }
    }
  }

  DC_recorder.close();
  spike_recorder.close();
  
  return 0;
} 
