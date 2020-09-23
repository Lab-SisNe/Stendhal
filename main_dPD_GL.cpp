/*
filename: main.cpp

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
#include <string>
#include <array>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <chrono>

// Include pcg if compiled with -D USE_PCG
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif

// Include stendhal headers
//#include "stendhalconfig.h"
#include "dPD_GL.hpp"
#include "gl_psc_exp.hpp"

/* Filename: main_dPD_GL.cpp
main function

*/

// Usage message
static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name << " seed t_sim" << std::endl
	    << "  seed: seed for random number generator (integer)" << std::endl
	    << "  t_sim: simulation time in (ms)" << std::endl
	    << "    if t_sim is ommitted, it defaults to 1000.0 ms"
	    << std::endl;
}

int main ( int argc, char* argv[] )
{
  double t_sim;
  if ((argc == 1) || (argc > 3)) {
    show_usage(argv[0]);
    return 1;
  }
  if (argc==2)
    t_sim = 1000.0;
  else
    t_sim = std::atof(argv[2]);
  int seed = std::atoi(argv[1]);
  std::cout << "seed: " << seed
	    << ", Simulation time: "
	    << t_sim << " ms" << std::endl;
  
  class stendhal::dPD_GL dpd_gl(seed);

  auto start = std::chrono::steady_clock::now();
  dpd_gl.prepare(); //call calibrate, create_pop and connect(void)
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end-start;
  std::cout << "prepare: " << diff.count() << " s\n";

  /*
  auto start = std::chrono::steady_clock::now();
  dpd_gl.calibrate();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> diff = end-start;
  std::cout << "calibrate: " << diff.count() << " s\n";

  start = std::chrono::steady_clock::now();
  dpd_gl.create_pop();
  end = std::chrono::steady_clock::now();
  diff = end-start;
  std::cout << "create_pop: " << diff.count() << " s" << std::endl;

  start = std::chrono::steady_clock::now();
  dpd_gl.connect("conn_NEST.txt");
  end = std::chrono::steady_clock::now();
  diff = end-start;
  std::cout << "connect: " << diff.count() << " s" << std::endl;
  */
  
  start = std::chrono::steady_clock::now();
  dpd_gl.simulate(t_sim);
  end = std::chrono::steady_clock::now();
  diff = end-start;
  std::cout << "simulate(" << t_sim << "): " << diff.count() << " s" << std::endl;
  
  return 0;
} 
