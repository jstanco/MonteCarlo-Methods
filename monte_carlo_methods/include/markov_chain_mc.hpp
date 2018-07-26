//created by john stanco on 5/10/18

// General purpose Markov-Chain Monte Carlo Simulation.
// Templated for any class with appropriate methods.
// Will recieve a state as input and evolve the state

#include "markov_chain.hpp"

#ifndef MARKOV_CHAIN_MC_H
#define MARKOV_CHAIN_MC_H


template<class T>
class markov_chain_mc : public markov_chain<T>{
public:
	markov_chain_mc(){}
	
	double expVal(const double (*f)(T*));
	double variance(const double (*f)(T*));

	void print_chain();
};


template<class T>
inline double
markov_chain_mc<T>::expVal(const double (*f)(T*)){
	double acc = 0;
	double len = this->chain.size();
	for(size_t i = 0; i < len; i++){
		acc += f(this->chain[i]);
	}
	return acc / len;
}


template<class T>
inline double
markov_chain_mc<T>::variance(const double (*f)(T*)){
	double acc = 0;
	double len = this->chain.size();
	double mean = this->expVal(f);
	for(size_t i = 0; i < len; i++){
		acc += pow(f(this->chain[i]) - mean, 2);
	}
	return acc / len;
}


template<class T>
inline void
markov_chain_mc<T>::print_chain(){
	for(size_t i = 0; i < this->chain.size(); i++){
		std::cout << i << "    " << this->chain[i].logProb() << std::endl;
	}
}


#endif /* MARKOV_CHAIN_MC_H */