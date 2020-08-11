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

#ifndef STENDHAL_PD_GL_HPP
#define STENDHAL_PD_GL_HPP

/*
Name: PD_GL - Implementation of the Potjans and Diesmann (1) microcircuit model
              using the Galves and Locherbach (2) neuron model.

References:
  (1) Tobias Potjans and Markus Diesmann (2014). 
  (2) Antonio Galve and Eva Locherbach (2013).
*/


namespace stendhal
{
  
  class PD_GL
  {
  public:
    PD_GL();
    ~PD_GL();

  private:
    // define simulation parameters
    struct sim_params_
    {
      double dt = 0.1;
      double t_sim = 1000.0;
      int n_threads = 1;
      int n_mpi = 1;
    };

    // define network params
    struct network_params_
    {
      // Number of neurons per layer
      // L23E, L23I, L4E, L4I, L5E, L5I, L6E, L6I
      double N_full[8] = {20683, 5834, 21915, 5479, 4850, 1065, 14395};
      // Connection probability
      double conn_prob[8][8] = {
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
      double K_ext[8] = {1600, 1500, 2100, 1900, 2000, 1200, 2900, 2100};
      // Factor to scale indegrees
      double K_scaling = 1.0;
      // Factor to scale the number of neurons
      double N_scaling = 1.0;
      // Mean postsynaptic current amplitude (in pA)
      double PSP_e = 87.8;
      // Relative inhibitory synaptic strength
      double g = -4.0;
      // Rate of Poisonian spike generator (in Hz)
      double bg_rage = 8.0;
      // Flag to turn ON or OFF poissonian background input
      // when false, use DC
      bool poisson_input = true;
      // mean delay of excitatory connections (in ms)
      double mean_delay_exc = 1.5;
      // mean delay of inhibitory connections (in ms)
      double mean_delay_inh = 0.8;
      // relative standard deviation of delays
      double rel_std_delay = 0.5;
    };
  
    // define neuron params
    struct neuron_params_ {
      // mean initial membrane potential (mV)
      double V0_mean = 7.0;
      // standard deviation of the initial membrane potential (V)
      double V0_sd = 10.0;
      // Reset membrane potential (mV)
      double E_L = 0.0;
      // LIF Threshold potential (mV)
      double V_th = 15.0;
      // Reset membrane potential (mV)
      double V_reset = 0.0;
      // Membrane capacitance (pF)
      double C_m = 250.0;
      // Membrane time constant (ms)
      double tau_m = 10.0;
      // Time constant of excitatory postsynaptic currents (ms)
      double tau_syn_ex = 0.5;
      // Time constant of inhibitory postsynaptic currents (ms)
      double tau_syn_in = 0.5;
      // Refractory preiod (ms)
      double tau_ref = 2.0;
      // Slope of the firing probability function (1/mV)
      double gamma = 0.1;
      // Curvature of the firing probability function (unitless)
      double r = 0.4;
      // Rheobase potential; potential in which firing probability > 0
      double V_rheo = 15.0;
    };
  
  public:
  }

} // namespace Standhal

#endif STENDHAL_PD_GL_HPP
