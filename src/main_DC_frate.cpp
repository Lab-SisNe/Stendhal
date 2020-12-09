/*
filename: main_DC_frate.cpp

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
	    << "  DC_start: start value of DC current:" << std::endl
	    << "  DC_end: final value of DC current:" << std::endl
	    << "  DC_step: step value of DC current:" << std::endl
	    << std::endl
	    << " Simulate N neurons to constant DC current input."
	    << std::endl;
}

int main ( int argc, char* argv[] )
{
  // Get input arguments
  if ((argc != 7)) {
    show_usage(argv[0]);
    return 1;
  }
  int seed = std::atoi(argv[1]);
  double t_sim = std::atof(argv[2]);
  double delta_t = std::atof(argv[3]);
  double DC_start = std::atof(argv[4]);
  double DC_end = std::atof(argv[5]);
  double DC_step = std::atof(argv[6]);
  std::cout << "Random Number Generator Engine: ";
  std::cout << "mt19937" << std::endl;
  std::cout << "seed: " << seed
	    << ", Simulation time: "
	    << t_sim << " ms"
	    << ", Simulation time step: "
	    << delta_t << " ms" << std::endl;
  std::cout << "DC_start = " << DC_start
	    << ", DC_end = " << DC_end
	    << ", DC_step = " << DC_step << std::endl;

  // neurons, DC, N
  std::vector<stendhal::gl_psc_exp*> neurons;
  std::vector<double> DC;
  size_t N {0};
  unsigned int buffer_pos {0};
  unsigned int buffer_size {5};
  double t = 0.0;
  
  // random number generators
  std::mt19937 rng;
  // uniform random number generator
  std::uniform_real_distribution<> udist;

  rng.seed(seed);
  
  // Timer
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::time_point<std::chrono::steady_clock> end;
  std::chrono::duration<double> diff;

  // outpuf file
  std::ofstream spike_recorder;
  std::string spike_recorder_file {"spike_recorder.txt"};

  // open output file
  spike_recorder.open(spike_recorder_file);
  spike_recorder << "# dt = " << delta_t << std::endl;
  spike_recorder << "# time,  neuron_ID,  V_m" << std::endl;
  
  // Create neurons
  for (double dc=DC_start; dc<=DC_end; dc+=DC_step) {
    N++; // increment neuron counter
    // append new gl_psc_exp instance
    neurons.push_back(new stendhal::gl_psc_exp(delta_t, &rng, &udist, &buffer_pos, &buffer_size));
    DC.push_back(dc);
  }

  // simulate
  double V_spiked;
  while (t<=t_sim) {
    t += delta_t;
    for (std::vector<stendhal::gl_psc_exp*>::iterator it=neurons.begin(); it!=neurons.end(); it++) {
      (*it)->add_DC_input(DC[(*it)->get_id()-1]);
      V_spiked = (*it)->evaluated();
      if (V_spiked > 0) {
	spike_recorder << (int)std::round(t/delta_t) << ", "
		       << (*it)->get_id() << ", "
		       << V_spiked << std::endl;
      }
    }
  }

  spike_recorder.close();
  
  return 0;
} 
