//created by John Stanco 7/27/18


/// User defined headers
#include "../include/ising.hpp"
#include "../include/graph.hpp"
#include "../include/markov_chain.hpp"


typedef Ising::wolff_ising<Graph::hypercubic_lattice_graph> square_ising;
typedef markov_chain_mc<square_ising> ising_mc;


void test_ising_mc(){
	square_ising init(1, 0, 1, {10, 10});

	ising_mc mc;
	timer t;



	t.start();
	square_ising rslt = mc.run(init, 1e3);
	t.stop();
	t.print();

	std::cout << mean(rslt.E_data()) << std::endl;
	std::cout << mean(rslt.M_data()) << std::endl;
}


int main(){

	test_ising_mc();

	return 1;
}