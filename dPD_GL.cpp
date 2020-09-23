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
// Include C++ standard headers
#include <fstream>
#include <cmath>

// Include stendhal headers
#include "dPD_GL.hpp"

namespace stendhal
{
  double dPD_GL::t = 0.0;
  //const int dPD_GL::N_layers = 8;
  
  // dPD_GL constructor
  dPD_GL::dPD_GL(int s)
  {
    // Set random number generator seed
    seed = s;
    
    struct simulation_parameters sim_params;
    struct network_parameters net_params;

    // Seed random number generator
    rng.seed(seed);

    // Calculate auxiliary parameters
    calibrate();
    
    // Initialize network
    // Create neuron population
    create_pop();
    
    // Create connection
    connect();

    // open output file;
    spike_recorder.open(spike_recorder_file);
    spike_recorder << "dt = " << sim_params.delta_t << '\n';
    spike_recorder << "time,  neuron_ID,  V_m\n";
    
  } // class dPD_DL constructor


  // dPD_GL destructor
  dPD_GL::~dPD_GL()
  {
    if (spike_recorder.is_open())
      spike_recorder.close();
  } // dPD_GL desctructor

  // Calculate auxiliary parameters
  void dPD_GL::calibrate(void)
  {
    // Calculate number of neurons based on scaling factor N_scaling
    N_cumsum[0]=0;
    for (int i=0; i<N_layers; i++) {
      // Scaled number of neurons
      N_scaled[i] = (unsigned int)std::round(net_params.N_full[i]*net_params.N_scaling);
      // Cumulative sum of number of neurons
      N_cumsum[i+1] = N_cumsum[i] + N_scaled[i];
      // population ID range per layer
      // note that neuron ID starts at 1 whereas
      // variable index starts at 0 (obvious)
      pop_ID[i][0] = N_cumsum[i]+1;
      pop_ID[i][1] = N_cumsum[i+1];
      // mean of the poisson input
      lam[i] = net_params.bg_rate * net_params.K_ext[i] * sim_params.delta_t * 1e-3;
      // DC input that would have the same effect as the poisson input
      DC[i] = net_params.bg_rate * net_params.K_ext[i] * net_params.PSP_e * net_params.tau_syn_ex * 1e-3;
    } // pop
    
    // Calculate number of synapses based on N_scaled and conn_prob
    // pre pop
    for (int i=0; i<N_layers; i++) {
      // post pop
      for (int j=0; j<N_layers; j++) {
	// calculate number of synapses
	double Ca = net_params.conn_prob[j][i];
	int Npre = N_scaled[i];
	int Npost = N_scaled[j];
	K_scaled[j][i] = (unsigned int)std::round(std::log(1.0-Ca)/std::log(1.0-(1.0/(double)(Npre*Npost))));
      } // post pop
    } // pre pop
  } // calibrate

  // create neuron population
  void dPD_GL::create_pop(void)
  {
    for (int i=0; i<N_layers; i++) {
      // Append neuron one by one to neurons vector
      // We must add one by one to correctly set node count and node ID
      for (int j=0; j<N_scaled[i]; j++) {
	neurons.push_back(new stendhal::gl_psc_exp(sim_params.delta_t, &rng, &udist, &buffer_pos, &buffer_size));
      }
    }
  } // create pop

  void dPD_GL::connect(void)
  {
    std::ofstream outfile ("connection.txt");
    std::uniform_int_distribution<>::param_type uintp;
    int pre_ID; // presynaptic ID
    int post_ID; // postsynaptic ID
    double w; // weight
    double d; // delay

    // save connection table to file
    // Format:
    //   N_layers:
    //   Neuron ID per layer: initial_ID final_ID
    // pre_ID post_ID weight delay
    if (!outfile.is_open())
      outfile.open("connection.txt", std::ofstream::out);
    outfile << "  N_layers: " << N_layers << '\n';
    for (int n=0; n<N_layers; n++)
      outfile << "  " << pop_ID[n][0] << ' ' << pop_ID[n][1] << '\n';
    outfile << '\n';
    
    // create connection
    // pre pop
    for (int i=0; i<N_layers; i++) {
      // post pop
      for (int j=0; j<N_layers; j++) {
	// number of synapses
	for (int n=0; n<K_scaled[j][i]; n++) {
	  // draw uniform integer between range pop_ID[i][0] and pop_ID[i][1]
	  uintp = std::uniform_int_distribution<>::param_type (pop_ID[i][0], pop_ID[i][1]);
	  pre_ID = udist_int(rng, uintp);
	  // draw uniform integer between range pop_ID[j][0] and pop_ID[j][1]
	  uintp = std::uniform_int_distribution<>::param_type(pop_ID[j][0], pop_ID[j][1]);
	  post_ID = udist_int(rng, uintp);
	  // draw normal distribution for synaptic weight and delay
	  if ((i%2) == 0) { // pre is excitatory
	    // weight
	    w = 0.0;
	    double w_ = net_params.PSP_e;
	    double w_sd = w_ * net_params.PSP_sd;
	    while (w <= 0) // excitatory weight must be positive
	      w = w_ + w_sd * ndist(rng);
	    // delay
	    d = 0.0;
	    double d_ = net_params.mean_delay_exc;
	    double d_sd = d_ * net_params.rel_std_delay;
	    while (d < sim_params.delta_t) // delay must be >= delta_t
	      d = d_ + d_sd * ndist(rng);
	  }
	  else { // pre is inhibitory
	    // weight
	    w = 0.0;
	    double w_ = net_params.PSP_e * net_params.g;
	    double w_sd = w_ * net_params.PSP_sd;
	    while (w >= 0) // inhibitory weight must be negative
	      w = w_ - w_sd * ndist(rng);
	    // delay
	    d = 0.0;
	    double d_ = net_params.mean_delay_inh;
	    double d_sd = d_ * net_params.rel_std_delay;
	    while (d < sim_params.delta_t) // delay must be >= delta_t
	      d = d_ + d_sd * ndist(rng);
	  }
	  neurons[pre_ID-1]->connect(neurons[post_ID-1], w, d);
	  outfile << pre_ID << ", "  << post_ID << ", " << w << ", " << d << '\n';
	} // loop for number of synapses
      } // loop for post-synaptic population
    } // loop for pre-synaptic population
    // close file
    outfile.close();
  } // connect

  // simulate
  void dPD_GL::simulate(double t_sim)
  {
    double V_spiked;
    // Reopen file in append mode if not open
    if (!spike_recorder.is_open())
      spike_recorder.open(spike_recorder_file, std::ofstream::out | std::ofstream::app);

    // Apply input
    unsigned int p;
    std::poisson_distribution<>::param_type pparam;
    // iterate through layers
    for (int n=0; n<N_layers; n++) {
      // set lambda for poisson distribution
      pparam = std::poisson_distribution<>::param_type (lam[n]);
      // iterate through pop_ID[n][0] to pop_ID[n][1]
      for (int i=pop_ID[n][0]; i<=pop_ID[n][1]; i++) {
	if (net_params.poisson_input) { // poisson input
	  p = pdist(rng, pparam);
	  neurons[i-1]->add_input(p * net_params.PSP_e, (unsigned int)0);
	}
	else { // DC input
	  neurons[i-1]->add_DC_input(DC[n]);
	}
      }
    }
    
    // evolve time
    t += sim_params.delta_t;

    // evaluate
    for (std::vector<gl_psc_exp*>::iterator it=neurons.begin(); it!=neurons.end(); it++) {
      V_spiked = (*it)->evaluated();  // returns V_m prior to spike, 0.0 otherwise
      // store spike output to file
      if (V_spiked > 0)
	spike_recorder << (int)std::round(t/sim_params.delta_t) << ", " << (*it)->get_id() << ", " << V_spiked << '\n';
    }

    // advance ring buffer position
    buffer_pos = (buffer_pos+1) % buffer_size;
    
  } // simulate
} // namespace stendhal

