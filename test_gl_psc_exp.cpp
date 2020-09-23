// Use option -D USE_PCG to define USE_PCG

// Includes from C++
#include <iostream>
#include <random>

// Include pcg if compiled with -D USE_PCG
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif

// Includs from stendhal
#include "gl_psc_exp.hpp"
#include "connection.hpp"


int main( int argc, char* argv[] )
{
#ifdef USE_PCG
  pcg32 rng;
#else
  std::mt19937 rng;
#endif

  rng.seed(1000);
  
  std::uniform_real_distribution<> udist;
  std::normal_distribution<double> ndist;
  std::poisson_distribution<> pdist;

  double delta_t = 0.1;

  unsigned int buff_size {50};
  unsigned int buff_pos {0};
  
  double bg_freq = 8.0 * 2000; // Hz
  double lam = bg_freq * delta_t * 1e-3;
  std::poisson_distribution<>::param_type pparam{lam};

  stendhal::gl_psc_exp neuron(delta_t, &rng, &udist, &buff_pos, &buff_size);

  std::vector<stendhal::gl_psc_exp *> neuron_v; // vector storing pointer to neurons

  for (int i=0; i<10; i++) {
    neuron_v.push_back(new stendhal::gl_psc_exp(delta_t, &rng, &udist, &buff_pos, &buff_size));
    std::normal_distribution<double>::param_type nparam{1500.0, 35.0};
    double w = ndist(rng,nparam);
    nparam = std::normal_distribution<double>::param_type{1.5, 0.75};
    double d = ndist(rng,nparam);
    neuron.connect(neuron_v[i], w, d);
  }

  std::vector<stendhal::gl_psc_exp> neuron_v2 (10, {delta_t, &rng, &udist, &buff_pos, &buff_size});

  for (std::vector<stendhal::gl_psc_exp>::iterator it=neuron_v2.begin(); it!=neuron_v2.end(); it++) {
    std::normal_distribution<double>::param_type nparam{1500.0, 35.0};
    double w = ndist(rng,nparam);
    nparam = std::normal_distribution<double>::param_type{1.5, 0.75};
    double d = ndist(rng,nparam);
    neuron.connect(&(*it), w, d);
  }
      
  double V_spiked;
  for (double t=0.0; t<100.0; t+=delta_t) {
    unsigned int poisson = pdist(rng, pparam);
    neuron.add_input(poisson*87.8,(unsigned int)0);
    V_spiked = neuron.evaluated();
    if (V_spiked > 0)
      std::cout << t << ", " << neuron.get_id() << ", " << V_spiked << "\n";
      //std::cout << neuron.get_id() << ": " << t << ", " << V_spiked << "\n";
    for (std::vector<stendhal::gl_psc_exp*>::iterator it=neuron_v.begin(); it!=neuron_v.end(); it++) {
      V_spiked = (*it)->evaluated();
      if (V_spiked > 0)
	std::cout << t << ", " << (*it)->get_id() << ", " << V_spiked << "\n";
    }

    //    std::cout << "test\n";
    for (std::vector<stendhal::gl_psc_exp>::iterator it=neuron_v2.begin(); it!=neuron_v2.end(); it++) {
      V_spiked = (*it).evaluated();
      if (V_spiked > 0)
	std::cout << t << ", " << (*it).get_id() << ", " << V_spiked << "\n";
	//std::cout << (*it).get_id() << ": " << t << ", " << V_spiked << "\n";
    }
    buff_pos = (++buff_pos) % buff_size;
  }
  std::cout << '\n';

  return 0;
}
