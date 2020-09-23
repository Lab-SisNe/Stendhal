#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#define USE_PCG
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif


// Define layer class
class layer_class
{
private:
  // define neuron params
  struct neuron_parameters
  {
    // mean initial membrane potential (mV)
    double V0_mean {7.0};
    // standard deviation of the initial membrane potential (V)
    double V0_sd {10.0};
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
    // Slope of the firing probability function (1/mV)
    double gamma {0.1};
    // Curvature of the firing probability function (unitless)
    double r {0.4};
    // Rheobase potential; potential in which firing probability > 0
    double V_rheo {15.0};
  }; // struct neuron_parameters

  // number of elements
  int size;
  // state variables
  double *V_m; // membrane potential
  double *I_syn_ex; // excitatory post-synaptic current
  double *I_syn_in; // inhibitory post-synaptic current
  double *I_ext; // external input current
  int *is_ref; // refractory counter

  // buffer
  std::vector<double> *W_ex;
  std::vector<double> *W_in;
      
public:
  layer_class(int N, int len);
  ~layer_class();

  void resize_buffer(int len)
  {
    for (int i=0; i<size; i++) {
      W_ex[i].resize(len, 0.0);
      W_in[i].resize(len, 0.0);
    }
  };

  void evaluate(double dt);
}; // class layer_class

// layer class constructor
layer_class::layer_class(int N, int len = 50):
  size(N)
{
  V_m = new double[N];
  I_syn_ex = new double[N];
  I_syn_in = new double[N];
  I_ext = new double[N];
  is_ref = new int[N];
  W_ex = new std::vector<double>[N];
  W_in = new std::vector<double>[N];
  for (int i=0; i<N; i++) {
    W_ex[i].resize(len, 0.0);
    W_in[i].resize(len, 0.0);
  }
}

layer_class::~layer_class()
{
  delete []V_m;
  delete []I_syn_ex;
  delete []I_syn_in;
  delete []I_ext;
  delete []is_ref;
  delete []W_ex;
  delete []W_in;
}

#ifdef USE_PCG
double draw(pcg32 *rng, std::poisson_distribution<> *dist)
#else
double draw(std::mt19937 *rng, std::poisson_distribution<> *dist)
#endif
{
  return (*dist)(*rng);
}
 
#ifdef USE_PCG
double draw(pcg32 *rng, std::uniform_real_distribution<> *dist)
#else
double draw(std::mt19937 *rng, std::uniform_real_distribution<> *dist)
#endif
{
  return (*dist)(*rng);
}

class pdraw
{
private:
#ifdef USE_PCG
  pcg32 *prng;
#else
  std::mt19937 *prng;
#endif
  std::uniform_real_distribution<> *pudist;

public:
  pdraw() :
    pdraw(0)
  {
  }
  
  pdraw(int seed_val)
  {
#ifdef USE_PCG
    prng = new pcg32;
#else
    prng = new std::mg19937;
#endif
    seed(seed_val);
    pudist = new std::uniform_real_distribution<>;
  }
#ifdef USE_PCG
  pdraw(pcg32 *rng, std::uniform_real_distribution<> *udist)
#else
  pdraw(std::mg19938 *rng, std::uniform_real_distribution<> *udist)
#endif
  {
    std::cout << "inside pdraw(*, *)\n";
    prng = rng;
    pudist = udist;
  }
  ~pdraw()
  {
    //    delete prng;
    //    delete pudist;
  }

  double draw()
  {
    return (*pudist)(*prng);
  }

  void seed(int seed)
  {
    prng->seed(seed);
  }

  void set(std::uniform_real_distribution<>::param_type p)
  {
    pudist->param(p);
  }

  double a(void)
  {
    return pudist->a();
  }

  double b(void)
  {
    return pudist->b();
  }

  std::uniform_real_distribution<>* dist_pointer(void)
  {
    return pudist;
  }
};

 
 
int main(int argc, char* argv[]) {
#ifdef USE_PCG
  pcg32 rng;
#else
  std::mt19937 rng;
#endif
  int seed {55};

  std::poisson_distribution<> rng_p;
  std::uniform_real_distribution<> udist;
  
  double dt = 0.1; // ms
  double freq = 8.0*1800; // Hz

  double lam = freq * dt * 1e-3;
  std::poisson_distribution<>::param_type p{lam};

  rng.seed(seed);
  for (int i=0; i<10; i++)
    std::cout << rng_p(rng, p) << ' ';
  std::cout << '\n';

  rng.seed(0);
  for (int i=0; i<10; i++)
    std::cout << rng_p(rng, p) << ' ';
  std::cout << '\n';  

  rng.seed(0);
  for (int i=0; i<10; i++)
    std::cout << udist(rng) << ' ';
  std::cout << '\n';  

  rng.seed(0);
  //  rng_p.param(p);
  for (int i=0; i<10; i++)
    std::cout << draw(&rng, &udist) << ' ';
  std::cout << '\n';

  pdraw pdraw1(&rng, &udist);
  for (int i=0; i<10; i++)
    std::cout << pdraw1.draw() << ' ';
  std::cout << '\n';

  pdraw pdraw2(0);
  //  pdraw2.set(udist.param());
  for (int i=0; i<10; i++)
    std::cout << pdraw2.draw() << ' ';
  std::cout << '\n';

  pdraw2.seed(1000);
  for (int i=0; i<10; i++)
    std::cout << pdraw2.draw() << ' ';
  std::cout << '\n';

  std::cout << "udist: " << udist.a() << ", " << udist.b() << '\n';
  std::cout << "pdraw1: " << pdraw1.a() << ", " << pdraw1.b() << '\n';
  std::cout << "pdraw2: " << pdraw2.a() << ", " << pdraw2.b() << '\n';

  pdraw pdraw3(0);
  std::cout << "pdraw3->pudist *: " << pdraw3.dist_pointer() << '\n';
  std::cout << "pdraw3: " << pdraw3.a() << ", " << pdraw3.b() << '\n';
  pdraw pdraw4;
  std::cout << "pdraw4->pudist *: " << pdraw4.dist_pointer() << '\n';
  std::cout << "pdraw4: " << pdraw4.a() << ", " << pdraw4.b() << '\n';
  
#ifdef USE_PCG
  pcg32 rng3;
#else
  std::mt19937 rng3;
#endif
  std::uniform_real_distribution<> udist3;
  std::cout << "udist3: " << udist3.a() << ", " << udist3.b() << '\n';

#ifdef USE_PCG
  pcg32 *prng4 = new pcg32;
#else
  std::mt19937 *prng4 = new std::mt19937;
#endif
  std::uniform_real_distribution<> *pudist4 = new std::uniform_real_distribution<>;
  std::cout << "pudist4: " << pudist4->a() << ", " << pudist4->b() << '\n';  
  
  //  unsigned char buff[2] {0X01,0X10};
  //  unsigned char buff[2] {1, 32};
  //  unsigned int buff {32};
  std::vector<bool> buff_vector;
  std::bitset<10> buff;
  //  std::cout << "buff[0]: " << buff[0] << '\n';
  //  std::cout << "buff[1]: " << buff[1] << '\n';
  std::cout << sizeof(buff) << '\n';

  
  return 0;
}
