//created by john stanco on 5/10/18

//General Purpose Metropolis Hastings Monte Carlo Simulation (Templated to recieve any type that contains basic methods

#include "markov_chain.hpp"

#ifndef metropolis_hastings_mc_h
#define metropolis_hastings_mc_h


template<class T>
class metropolis_hastings_mc : public markov_chain<T>{
public:
	metropolis_hastings_mc(){}
	
	double expVal(const double (*f)(T*));
	double variance(const double (*f)(T*));

	void printChain();
};


template<class T>
inline double
metropolis_hastings_mc<T>::expVal(const double (*f)(T*))
{
	double acc = 0;
	double len = this->chain.size();
	for(unsigned i = 0; i < len; i++)
	{
		acc += f(this->chain[i]);
	}
	return acc / len;
}


template<class T>
inline double
metropolis_hastings_mc<T>::variance(const double (*f)(T*)){
	double acc = 0;
	double len = this->chain.size();
	double mean = this->expVal(f);
	for(unsigned i = 0; i < len; i++)
	{
		acc += pow(f(this->chain[i]) - mean, 2);
	}
	return acc / len;
}


template<class T>
inline void
metropolis_hastings_mc<T>::printChain(){
	for(unsigned i = 0; i < this->chain.size(); i++)
	{
		std::cout << i << "    " << this->chain[i].logProb() << std::endl;
	}
}


#endif /* metropolis_hastings_mc_h */