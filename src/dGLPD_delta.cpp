/*
filename: dGLPD_delta.cpp

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
// Include C++ standard headers
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
//#include <chrono>

// Include stendhal headers
#include "dGLPD_delta.hpp"

namespace stendhal
{
  double dGLPD::t = 0.0;
  int dGLPD::nt = 0;
  //const int dGLPD::N_layers = 8;
  
  // dGLPD constructor
  dGLPD::dGLPD(int s, double delta_t, bool arec)
  {
    // Set random number generator seed
    seed = s;
    
    // Seed random number generator; defaults to 55
    prng.seed(seed);

    // Set simlation time step; defaults to 0.1 ms
    sim_params.delta_t = delta_t;
    // set ticks_per_ms based on delta_t value;
    sim_params.ticks_per_ms = std::round(1/delta_t);

    // open output file;
    spike_recorder.open(spike_recorder_file);
    spike_recorder << "# seed = " << seed << std::endl;
    spike_recorder << "# dt = " << sim_params.delta_t << std::endl;    
    spike_recorder << "# neuron_ID, step,  V_m (mV)" << std::endl;
    if (arec) {
      analog_rec = arec;
      analog_recorder.open(analog_recorder_file);
      analog_recorder << "# seed = " << seed << std::endl;
      analog_recorder << "# dt = " << sim_params.delta_t << std::endl;
      analog_recorder << "# neuron_ID, step, V_m (mV), ePSC (pA), iPSC (pA), I_ext (pA), ePSP (mV), iPSP (mV), V_ext (mV)" << std::endl;
    }
      
    
  } // class dPD_DL constructor

  // dGLPD destructor
  dGLPD::~dGLPD()
  {
    if (spike_recorder.is_open())
      spike_recorder.close();
    if (analog_recorder.is_open())
      analog_recorder.close();
  } // dGLPD desctructor

  // Calculate auxiliary parameters
  void dGLPD::calibrate(void)
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
      DC[i] = net_params.bg_rate * net_params.K_ext[i] * net_params.PSP_e * net_params.C_m * 1e-3;
    } // pop
    
    // Calculate number of synapses based on N_scaled and conn_prob
    // post pop
    for (int i=0; i<N_layers; i++) {
      // pre pop
      for (int j=0; j<N_layers; j++) {
	// calculate number of synapses
	double Ca = net_params.conn_prob[i][j];
	double Npre = (double)N_scaled[j];
	double Npost = (double)N_scaled[i];
	K_scaled[i][j] = (unsigned int)std::round(std::log(1.0-Ca)/std::log(1.0-1.0/(Npre*Npost)));
	// calculate weight matrix
	if ((j%2) == 0) {
	  weight_matrix[i][j] = net_params.PSP_e;
	  if ((j==2) && (i==0)) {
	    weight_matrix[i][j] = 2 * net_params.PSP_e;
	  }
	  delay_matrix[i][j] = net_params.mean_delay_exc;
	}
	else {
	  weight_matrix[i][j] = net_params.g * net_params.PSP_e;
	  delay_matrix[i][j] = net_params.mean_delay_inh;
	}
      } // pre pop
    } // post pop
  } // calibrate

  // create neuron population
  void dGLPD::create_pop(void)
  {
    for (int i=0; i<N_layers; i++) {
      // Append neuron one by one to neurons vector
      // We must add one by one to correctly set node count and node ID
      for (int j=0; j<(int)N_scaled[i]; j++) {
	neurons.push_back(new stendhal::gl_psc_delta(sim_params.delta_t, &prng, &udist, &buffer_pos, &buffer_size));
      }
    }
  } // create pop

  // Create connection
  void dGLPD::connect(void)
  {
    std::ofstream outfile ("connection.txt");
    std::uniform_int_distribution<>::param_type uintpre, uintpost;
    int pre_ID; // presynaptic ID
    int post_ID; // postsynaptic ID
    double w, w_, w_sd; // weight
    double d, d_, d_sd; // delay
    double d_max = 0.0;
    
    // save connection table to file
    // Format:
    //   N_layers:
    //   Neuron ID per layer: initial_ID final_ID
    // pre_ID post_ID weight delay
    if (!outfile.is_open())
      outfile.open("connection.txt", std::ofstream::out);
    outfile << "# N_layers: " << N_layers << '\n';
    for (int n=0; n<N_layers; n++)
      outfile << "# " << pop_ID[n][0] << ' ' << pop_ID[n][1] << '\n';
    outfile << '\n';
    
    // create connection
    // pre pop
    for (int j=0; j<N_layers; j++) {
      // post pop
      for (int i=0; i<N_layers; i++) {
	// uniform integer distribution parameters
	uintpre = std::uniform_int_distribution<>::param_type (pop_ID[j][0], pop_ID[j][1]); // post-synaptic parameters
	uintpost = std::uniform_int_distribution<>::param_type (pop_ID[i][0], pop_ID[i][1]); // post-synaptic parameters
	// weight parameters
	w_ = weight_matrix[i][j];
	w_sd = w_ * net_params.PSP_sd;
	// delay parameters
	d_ = delay_matrix[i][j];
	d_sd = d_ * net_params.rel_std_delay;
	// number of synapses
	for (int n=0; n<(int)K_scaled[i][j]; n++) {
	  // draw uniform integer between range pop_ID[j][0] and pop_ID[j][1]
	  pre_ID = udist_int(prng, uintpre);
	  // draw uniform integer between range pop_ID[i][0] and pop_ID[i][1]
	  post_ID = udist_int(prng, uintpost);
	  // draw normal distribution for synaptic weight and delay
	  // weight
	  w = 0.0;
	  if ((j%2) == 0) { // pre is excitatory
	    while (w <= 0) // excitatory weight must be positive
	      w = w_ + w_sd * ndist(prng);
	  }
	  else { // pre is inhibitory
	    while (w >= 0) // inhibitory weight must be negative
	      w = w_ - w_sd * ndist(prng);
	  }
	  // delay
	  d = 0.0;
	  while (d < sim_params.delta_t) // delay must be >= delta_t
	    d = d_ + d_sd * ndist(prng);

	  // update d_max to define buffer size
	  if (d > d_max)
	    d_max = d;
	  // create connetion
	  neurons[pre_ID-1]->connect(neurons[post_ID-1], w, d, std::round(d*sim_params.ticks_per_ms));
	  // store connection
	  outfile << pre_ID << ", "  << post_ID << ", " << w << ", " << d << '\n';
	} // loop for number of synapses
      } // loop for post-synaptic population
    } // loop for pre-synaptic population
    // close file
    outfile.close();
    // update buffer length if necessary
    if (d_max > (buffer_size*sim_params.delta_t))
      update_buffer_size(d_max);
  } // connect

  // Create connection from file
  void dGLPD::connect(std::string fname)
  {
    std::ifstream input_file(fname);
    double d_max = 0.0;

    if(!input_file.is_open())
      throw std::runtime_error("Could not open connection file");

    std::string line, word;
    std::vector<std::string> row;
    int i, j;
    double w, d;

    // read lines
    while (std::getline(input_file, line)) {
      if (line[0] == ' ' || line[0] == '#' || line.length() == 0) {
	continue;
      }
      // create string stream from current line
      std::stringstream ss(line);
      // clear row vector
      row.clear();
      // read each comma separated val in row vector
      while (std::getline(ss, word, ','))
	row.push_back(word);
      // convert each word to respective variables; 
      i = std::stoi(row[0]);
      j = std::stoi(row[1]);
      w = std::stof(row[2]);
      d = std::stof(row[3]);
      // Connect node
      neurons[i-1]->connect(neurons[j-1], w, d, std::round(d*sim_params.ticks_per_ms));
      if (d > d_max)
	d_max = d;
    }
    // update buffer length if necessary
    if (d_max > ((buffer_size)*sim_params.delta_t))
      update_buffer_size(d_max);
  } // connect from file


  // Prepare method
  // call: calibrate, create_node and connect
  void dGLPD::prepare()
  {
    // Calculate auxiliary parameters
    calibrate();
    
    // Initialize network
    // Create neuron population
    create_pop();
    
    // Create connection
    connect();
  } // prepare
  
  // simulate
  void dGLPD::simulate(double t_sim)
  {
    simulate_(std::round(t_sim*sim_params.ticks_per_ms));
  }

  void dGLPD::simulate_(int n_sim)
  {
    double V_spiked;
    unsigned int p;
    std::poisson_distribution<>::param_type pparam;

    // Reopen file in append mode if not open
    if (!spike_recorder.is_open())
      spike_recorder.open(spike_recorder_file, std::ofstream::out | std::ofstream::app);

    // iterate for a period of t_sim; to do so, t must be added to t_sim
    // to account for simulation starting at time t
    n_sim += nt;
    while (nt<n_sim) {
      // evolve time
      t += sim_params.delta_t; // ms
      nt++; // step

      // check if record analog data
      bool is_analog_rec = (analog_rec and ((nt % sim_params.ticks_per_ms) == 0) );

      // iterate through layers
      for (int n=0; n<N_layers; n++) {
	// set lambda for poisson distribution; shared between neurons of the same layer
	pparam = std::poisson_distribution<>::param_type (lam[n]);
	// iterate through pop_ID[n][0] to pop_ID[n][1]
	for (int i=pop_ID[n][0]; i<=(int)pop_ID[n][1]; i++) {
	  // Apply input
	  if (net_params.poisson_input) { // poisson input
	    p = pdist(prng, pparam); // draw poisson distribution with rate lam[n]
	    neurons[i-1]->add_input(p * net_params.PSP_e, (unsigned int)0);
	  }
	  else { // DC input
	    neurons[i-1]->add_DC_input(DC[n]);
	  } // external input
    
	  // evaluate
	  V_spiked = neurons[i-1]->evaluate();  // returns V_m when the neuron spiked; 0.0 otherwise
	  // store spike output to file when neuron spiked
	  if (V_spiked != 0)
	    spike_recorder << neurons[i-1]->get_id() << ", " << t << ", " << V_spiked << '\n';
	  // store analog data (membrane potentials and currents) to file
	  if (is_analog_rec) {
	    analog_recorder << neurons[i-1]->get_id() // neuron ID
			    << ", " << t // time step
			    << ", " << neurons[i-1]->get_Vm() // membrane potential
			    << ", " << neurons[i-1]->get_ePSC() // excitatory post synaptic current (ePSC) (mV)
			    << ", " << neurons[i-1]->get_iPSC() // inhibitory post synaptic current (iPSC) (mV)
			    << ", " << neurons[i-1]->get_I_ext() // external current (pA)
			    << ", " << neurons[i-1]->get_ePSP() // excitatory post synaptic potential (ePSP) (mV)
			    << ", " << neurons[i-1]->get_iPSP() // inhibitory post synaptic potential (iPSP) (mV)
			    << ", " << neurons[i-1]->get_V_ext() // external current induced potential change (mV)
			    << std::endl;
	  } // if is_analog_rec
	} // population
      } // layer

      // advance ring buffer position
      // all nodes share this same adress
      buffer_pos = (buffer_pos+1) % buffer_size;
    } // while loop    
  } // simulate

  // update buffer size of all nodes
  void dGLPD::update_buffer_size(double d)
  {
    // update buffer size according to d
    // set buffer_size to number os steps necessary to fit delay (d)
    // plus one, as during simulation the position buff_pos is at t+dt
    // if a delay of buff_size is to be added, it will be added at the
    // same position as t+dt, which will be erased after neuron evaluation.
    buffer_size = std::round(d*sim_params.ticks_per_ms)+1;
    // iterate through all neurons to update buffer size
    for (std::vector<gl_psc_delta*>::iterator it=neurons.begin(); it!=neurons.end(); it++)
      (*it)->resize_buffer(buffer_size);
  } // update buffer size

  // temporary
  double dGLPD::get_conn_prob(int i, int j)
  {
    return net_params.conn_prob[j][i];
  }

  int dGLPD::get_N_layers(void)
  {
    return N_layers;
  }

} // namespace stendhal

