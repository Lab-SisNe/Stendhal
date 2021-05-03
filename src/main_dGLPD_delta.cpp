/*
filename: main_dGLPD_delta.cpp

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

Aug 11, 2020
*/

// Include C++ standard headers
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <chrono>

// Include stendhal headers
#include "stendhalconfig.h"
#include "dGLPD_delta.hpp"
#include "gl_psc_delta.hpp"

/* Filename: main_dGLPD_delta.cpp
main function

*/

// Usage message
static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name << " seed t_sim <conn_csv_data>" << std::endl
	    << "  seed: seed for random number generator (integer)" << std::endl
	    << "  t_sim: simulation time in (ms)" << std::endl
	    << "  delta_t: simulation time step (ms)" << std::endl
	    << "  conn_csv_data (optional): CSV file containing connectivity data:" << std::endl
	    << "    if conn_csv_data is ommited, connection is randomly generated"
	    << std::endl;
}

int main ( int argc, char* argv[] )
{
  if ((argc < 4) || (argc > 5)) {
    show_usage(argv[0]);
    return 1;
  }
  int seed = std::atoi(argv[1]);
  double t_sim = std::atof(argv[2]);
  double delta_t = std::atof(argv[3]);
  int n_sim = std::round(t_sim/delta_t);
  int ticks_per_ms = std::round(1/delta_t);
  
  std::cout << "Random Number Generator Engine: ";
  std::cout << "xoroshiro128+" << std::endl;
  std::cout << "seed: " << seed
	    << ", Simulation time: "
	    << t_sim << " ms"
	    << ", Simulation time step: "
	    << delta_t << " ms";
  if (argc==5)
    std::cout << ", CSV_file: " << argv[4];
  std::cout << std::endl;

  class stendhal::dGLPD dglpd(seed, delta_t, false);

  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::time_point<std::chrono::steady_clock> end;
  std::chrono::duration<double> diff;
  
  if (argc != 5) {
    start = std::chrono::steady_clock::now();
    dglpd.prepare(); //call calibrate, create_pop and connect(void)
    end = std::chrono::steady_clock::now();
    diff = end-start;
    std::cout << "prepare: " << diff.count() << " s\n";
  }
  else {
    start = std::chrono::steady_clock::now();
    dglpd.calibrate();
    end = std::chrono::steady_clock::now();
    diff = end-start;
    std::cout << "calibrate: " << diff.count() << " s\n";

    start = std::chrono::steady_clock::now();
    dglpd.create_pop();
    end = std::chrono::steady_clock::now();
    diff = end-start;
    std::cout << "create_pop: " << diff.count() << " s" << std::endl;

    start = std::chrono::steady_clock::now();
    dglpd.connect(argv[4]);
    end = std::chrono::steady_clock::now();
    diff = end-start;
    std::cout << "connect: " << diff.count() << " s" << std::endl;
  }

  std::cout << "t_sim: " << t_sim << std::endl;
  start = std::chrono::steady_clock::now();
  dglpd.simulate(t_sim);
  end = std::chrono::steady_clock::now();
  diff = end-start;
  std::cout << "simulate(" << t_sim << "): " << diff.count() << " s" << std::endl;

  return 0;
} 
