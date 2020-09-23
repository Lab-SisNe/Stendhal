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

// Include stendhal headers
#include "stendhalconfig.h"
//#include "PD_GL.hpp"
/* Filename: main.cpp
main function

*/

// uncomment to use pcg random number generator
#define USE_PCG
//#ifdef USE_PCG
#include "pcg_cpp/pcg_random.hpp"
#endif
#include <random>

// Connector structure; stores target, weight and delay
struct connector
{
  int target;
  double weight;
  double delay;
};

// Discrete Galves and Locherbach model
class dGL
{
  // paramter key enumerator for use in swich case
  enum param_key {
		  invalid=-1,
		  delta_t,
		  E_L,
		  V_th,
		  V_reset,
		  C_m,
		  tau_m,
		  tau_syn_ex,
		  tau_syn_in,
		  tau_ref,
		  I_0,
		  gamma,
		  r,
		  V_rheo
  };

  // return enumerator from string
  param_key haskey (std::string const& inString)
  {
    if (inString == "delta_t") return delta_t;
    if (inString == "E_L") return E_L;
    if (inString == "V_th") return V_th;
    if (inString == "V_reset") return V_reset;
    if (inString == "C_m") return C_m;
    if (inString == "tau_m") return tau_m;
    if (inString == "tau_syn_ex") return tau_syn_ex;
    if (inString == "tau_syn_in") return tau_syn_in;
    if (inString == "tau_ref") return tau_ref;
    if (inString == "I_0") return I_0;
    if (inString == "gamma") return gamma;
    if (inString == "r") return r;
    if (inString == "V_rheo") return V_rheo;
    return invalid;
  }

  // random number generator
  std::mt19937 gen;
  std::uniform_real_distribution<> urand{0.0, 1.0};
  
  /*
  // Dictionary storing neuron parameters
  std::map <std::string, double> neuron_params {
    // Reset membrane potential (mV)
    {"E_L", 0.0},
    // LIF Threshold potential (mV)
    // {"V_th", 15.0},
    // Reset membrane potential (mV)
    {"V_reset", 0.0},
    // Membrane capacitance (pF)
    {"C_m", 250.0},
    // Membrane time constant (ms)
    {"tau_m", 10.0},
    // Excitatory synaptic time constant (ms)
    {"tau_syn_ex", 0.5},
    // Inhibitory synaptic time constant (ms)
    {"tau_syn_in", 0.5},
    // Refractory period (ms)
    {"tau_ref", 2.0},
    // Slope of firing probability (1/mV)
    {"gamma", 0.1},
    // Curvature of firing probability (unitless)
    {"r", 0.4},
    // Rehobase potential, potential where firing probability becomes > 0
    {"V_rheo", 15.0}
  };
  */

  // Neuron paramters
  struct neuron_parameters
  {
    // Simulation time step (ms)
    double delta_t {0.1};
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
    // Constant input current
    double I_0 {0.0};
    // Slope of the firing probability function (1/mV)
    double gamma {0.1};
    // Curvature of the firing probability function (unitless)
    double r {0.4};
    // Rheobase potential; potential in which firing probability > 0
    double V_rheo {15.0};
  } params;

  // State variables
  struct state_var
  {
    // Membrane potential (mV)
    double V_m {0.0};
    // Excitatory synaptic current (pA)
    double I_syn_ex {0.0};
    // Inhibitory synaptic current (pA)
    double I_syn_in {0.0};
    // External current (pA)
    double I_ext {0.0};
    // Refractory step counter;
    int is_ref {0};
  } state;

  // Ring input buffer
  struct buffer
  {
    // buffer size
    int N;
    // current buffer position
    int curr_pos;
    // Excitatory input buffer (pA)
    std::vector<double> W_ex;
    // Inhibitory input buffer (pA)
    std::vector<double> W_in;
    // buffer initializer (N=50)
    buffer():
      N(50),
      curr_pos(0)
    {
      W_ex.resize(N,0.0);
      W_in.resize(N,0.0);
    }
  } buff;

  // Auxiliary coefficient
  struct auxiliary_coefficient
  {
    // membrane potential leakage factor (unitless)
    double rho_m;
    // excitatory post-synaptic current decay factor (unitless)
    double rho_ex;
    // inhibitory post-synaptic current decay factor (unitless)
    double rho_in;
    // excitatory current propagator (Ohm; resistance)
    double Gamma_ex;
    // inhibitory current propagator (Ohm; resistance)
    double Gamma_in;
    // External step current input propagator (Ohm; resistance)
    double Gamma_DC;
    // Refractory period duration (in steps)
    int ref_count;
  } aux;

  // Calculate auxiliary coefficients
  // Must be called whenever some of the parameter are updated,
  // or at initialization
  void calc_aux_coeff(void)
  {
    aux.rho_m = exp( -params.delta_t/params.tau_m );
    aux.rho_ex = exp( -params.delta_t/params.tau_syn_ex );
    aux.rho_in = exp( -params.delta_t/params.tau_syn_in );
    if (params.tau_m != params.tau_syn_ex)
      aux.Gamma_ex = \
	(aux.rho_m - aux.rho_ex) * params.tau_m * params.tau_syn_ex \
	/ ( params.tau_m - params.tau_syn_ex) / params.C_m;
    else
      aux.Gamma_ex = params.delta_t * aux.rho_m / params.C_m;
    if (params.tau_m != params.tau_syn_in)
      aux.Gamma_in = \
	(aux.rho_m - aux.rho_in) * params.tau_m * params.tau_syn_in \
	/ ( params.tau_m - params.tau_syn_in) / params.C_m;
    else
      aux.Gamma_in = params.delta_t * aux.rho_m / params.C_m;
    aux.Gamma_DC = params.tau_m * (1 - aux.rho_m) / params.C_m;
    aux.ref_count = (int)params.tau_ref / params.delta_t;
  }

public:
  dGL(int seed)
  {
    // calculate auxiliary coefficients
    calc_aux_coeff();
    //    resize_buff(50);
    // seed the random number generator with seed
    gen.seed(seed);
  }

  ~dGL()
  {
  }
  
  // return buffer size
  int get_buff_size(void)
  {
    return buff.N;
  }

  // resize buffer
  void resize_buff(int N)
  {
    if (N != buff.N) {
      buff.W_ex.resize(N, 0.0);
      buff.W_in.resize(N, 0.0);
    }
  }
  
  // resize buffer from max_delay in (ms)
  void resize_buff(double max_delay)
  {
    resize_buff( (int)(max_delay/params.delta_t) );
  }

  // get current W_ex value
  double get_w_ex(void)
  {
    return buff.W_ex[buff.curr_pos];
  }

  // get current W_in value
  double get_w_in(void)
  {
    return buff.W_in[buff.curr_pos];
  }

  // Advance buffer position
  void advance_buffer(void)
  {
    buff.curr_pos++;
    if (buff.curr_pos >= buff.N)
      buff.curr_pos-=buff.N;
  }

  // set parameters
  void set_params(std::string key, double value)
  {
    /*
    try {
      neuron_params.at(key) = value;
    }
    catch (const std::out_of_range& e) {
      std::cout << "Invalid key: " << key << "\n";
    }
    */
    switch (haskey(key)) {
    case delta_t:
      params.delta_t = value;
      break;
    case E_L:
      params.E_L = value;
      break;
    case V_th:
      params.V_th = value;
      break;
    case V_reset:
      params.V_reset = value;
      break;
    case C_m:
      params.C_m = value;
      break;
    case tau_m:
      params.tau_m = value;
      break;
    case tau_syn_ex:
      params.tau_syn_ex = value;
      break;
    case tau_syn_in:
      params.tau_syn_in = value;
      break;
    case tau_ref:
      params.tau_ref = value;
      break;
    case I_0:
      params.I_0 = value;
      break;
    case gamma:
      params.gamma = value;
      break;
    case r:
      params.r = value;
      break;
    case V_rheo:
      params.V_rheo = value;
      break;
    default:
      std::cout << "Invalid key (set_params): " << key << "\n";
      break;
    }
    // recalculate auxiliary coefficients when necessary
    if ((key == "tau_m") || (key == "tau_syn_ex") || (key == "tau_syn_in")) {
      calc_aux_coeff();
    }
    if (key == "C_m") {
      if (params.tau_m != params.tau_syn_ex)
	aux.Gamma_ex = \
	  (aux.rho_m - aux.rho_ex) * params.tau_m * params.tau_syn_ex \
	  / ( params.tau_m - params.tau_syn_ex) / params.C_m;
      else
	aux.Gamma_ex = params.delta_t * aux.rho_m / params.C_m;
      if (params.tau_m != params.tau_syn_in)
	aux.Gamma_in = \
	  (aux.rho_m - aux.rho_in) * params.tau_m * params.tau_syn_in \
	  / ( params.tau_m - params.tau_syn_in) / params.C_m;
      else
	aux.Gamma_in = params.delta_t * aux.rho_m / params.C_m;
      aux.Gamma_DC = params.tau_m * (1 - aux.rho_m) / params.C_m;
    }
  }
  // set parameters from std::pair
  void set_params(std::pair <std::string, double> *par_pair)
  {
    set_params(par_pair->first, par_pair->second);
  }
  // set parameters from std::map
  void set_params(std::map <std::string, double> *par_map)
  {
    for (std::map<std::string, double>::iterator it = par_map->begin(); it != par_map->end(); ++it)
      set_params(it->first, it->second);
  }

  double get_params(std::string key)
  {
    switch (haskey(key)) {
    case delta_t:
      return params.delta_t;
    case E_L:

      return params.E_L;
    case V_th:
      return params.V_th;
    case V_reset:
      return params.V_reset;
    case C_m:
      return params.C_m;
    case tau_m:
      return params.tau_m;
    case tau_syn_ex:
      return params.tau_syn_ex;
    case tau_syn_in:
      return params.tau_syn_in;
    case tau_ref:
      return params.tau_ref;
    case I_0:
      return params.I_0;
    case gamma:
      return params.gamma;
    case r:
      return params.r;
    case V_rheo:
      return params.V_rheo;
    default:
      std::cout << "Invalid key (get_params): " << key << "\n";
      return 0.0;
    }
  }

  // Calculate firing probability
  double phi(void)
  {
    double V_diff = state.V_m - params.V_rheo;
    if (V_diff < 0)
      V_diff = 0.0;
    return std::pow(params.gamma * V_diff, params.r);
  }
  
  // Evolve neuron dynamics to next time step;
  // returns true if neuron fired, false otherwise
  bool evaluate(void)
  {
    // return value; true if neuron fired, false otherwise
    bool ret {false};
    
    // update membrane potential if neuron not refractory
    if (aux.ref_count <= 0) {
      state.V_m = aux.rho_m * state.V_m + aux.Gamma_ex * state.I_syn_ex + \
	aux.Gamma_in * state.I_syn_in + aux.Gamma_DC * (state.I_ext + params.I_0);
      
      // read synaptic inputs and advance buffer position
      //    double w_ex = get_w_ex();
      //    double w_in = get_w_in();
      //    advance_buffer();
      // update synaptic current inputs
      state.I_syn_ex = aux.rho_ex * state.I_syn_ex + buff.W_ex[buff.curr_pos];
      state.I_syn_in = aux.rho_in * state.I_syn_in + buff.W_in[buff.curr_pos];
      
      // advance buffer position
      buff.curr_pos++;
      if (buff.curr_pos >= buff.N)
	buff.curr_pos-=buff.N;

      // Calculate firing probability
      double p = phi();
      // draw a uniform random number
      double u = urand(gen);
      // if p >= u, neuron fires
      if (p >= u) {
	// Set return value to true
	ret = true;
	// reset membrane potential
	state.V_m = params.V_reset;
	// reset input currents;
	state.I_syn_ex = 0.0;
	state.I_syn_in = 0.0;
	// set refractory period counter
	state.is_ref = aux.ref_count;
      }
    }
    else
      // decrease refractory counter
      aux.ref_count--;
    // return firing state
    return ret;
  }
};

// Potjans and Diesmann network  
class dPD
{
private:
  struct params {
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
    double conn_prob[8][8] {
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
    double K_ext[8] {1600, 1500, 2100, 1900, 2000, 1200, 2900, 2100};
    // Factor to scale indegrees
    double K_scaling {1.0};
    // Factor to scale the number of neurons
    double N_scaling {0.1};
    // Mean postsynaptic potential amplitude (in mV)
    double PSP_e {0.15};
    // Relative inhibitory synaptic strength
    double g {-4.0};
    // Rate of Poisonian spike generator (in Hz)
    double bg_rage {8.0};
    // Flag to turn ON or OFF poissonian background input
    // when false, use DC
    bool poisson_input {true};
    // mean delay of excitatory connections (in ms)
    double mean_delay_exc {1.5};
    // mean delay of inhibitory connections (in ms)
    double mean_delay_inh {0.8};
    // relative standard deviation of delays
    double rel_std_delay {0.5};
  };

public:
  dPD(void);
  ~dPD(void);
  
};

class layer
{
  //private:
public:
  int size;
  double *V_m;
  double *I_syn_ex;
  double *I_syn_in;
  double *I_ext;
  int *is_ref;

  std::vector<double> *W_ex;
  std::vector<double> *W_in;
  //public:
  layer(int N, int len):
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
      std::cout << "resize\n";
      W_ex[i].resize(len);
      W_in[i].resize(len);
    }
  };
  ~layer()
  {
    delete []V_m;
    delete []I_syn_ex;
    delete []I_syn_in;
    delete []I_ext;
    delete []is_ref;
    delete []W_ex;
    delete []W_in;
  };
};

int main ( int argc, char* argv[] )
{
  // Create the network
  layer l1(10,50);
  layer l2(30,50);

  l1.V_m[5]=10.0;
  for (int i=0; i<l1.size; i++)
    std::cout << l1.V_m[i] << ' ';
  std::cout << '\n';
  for (int i=0; i<l2.size; i++)
    std::cout << l2.V_m[i] << ' ';
  std::cout << '\n';
  std::cout << "W_ex size: " << l2.W_ex[0].size() << '\n';
  //  gl_neuron.set_params("V_th", 0.0);
  //  std::cout << gl_neuron.get_params("V_th") << '\n';
  
  return 0;
}
