// Use option -D USE_PCG to define USE_PCG

#include <iostream>
#ifdef USE_PCG
#include "pcg-cpp/pcg_random.hpp"
#endif
#include <random>

#include "simulation-manager.hpp"
#include "dummy.hpp"
#include "connection.hpp"

int main( int argc, char* argv[] )
{
  stendhal::SimulationManager sim_man;

  std::cout << "t: " << sim_man.get_time() << ", ";
  std::cout << "N_nodes: " << sim_man.get_n_nodes() << '\n';

  for (int n=0; n<5; n++)
    sim_man.append_node(new stendhal::dummy);

  std::cout << "t: " << sim_man.get_time() << ", ";
  std::cout << "N_nodes: " << sim_man.get_n_nodes() << '\n';
  
  for (int i=0; i<10; i++) {
    std::cout << "t: " << sim_man.get_time() << ", ";
    std::cout << "N_nodes: " << sim_man.get_n_nodes() << '\n';
    for (int n=0; n<sim_man.get_n_nodes(); n++) {
      std::cout << "    ";
      std::cout << sim_man.get_node(n+1)->get_node_ID();
      std::cout << ", ";
      std::cout << sim_man.get_node(n+1)->get_time();
      std::cout << '\n';
    }
  }

  return 0;
}
