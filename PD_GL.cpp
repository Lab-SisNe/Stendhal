/*
filename: PD_GL.cpp

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

#include <iostream>
#include "dPD_GL.hpp"

namespace stendhal
{

  PD_GL_::PD_GL_()
  {
    struct simulation_parameters sim_params;
    struct network_parameters net_params;
    
  }

  PD_GL_::~PD_GL_()
  {
  }
  
  void PD_GL_::print_N_full(void)
  {
    for (int n: net_params.N_full)
      std::cout << n << '\n';
  }
    
  void PD_GL_::print_conn_prob(int i, int j)
  {
    std::cout << i << ", " << j << ": " << net_params.conn_prob[i][j] << std::endl;
  }

} // namespace stendhal

