//created by john stanco on 5/10/18

//General Purpose Metropolis Hastings Monte Carlo Simulation (Templated to recieve any type that contains basic methods

#include "MarkovChain.hpp"

#ifndef Metropolis_Hastings_MC_h
#define Metropolis_Hastings_MC_h


template<class T>
class Metropolis_Hastings_MC : public MarkovChain<T>{
public:
	Metropolis_Hastings_MC(){}
	
	double expVal(const double (*f)(T*));
	double variance(const double (*f)(T*));

	void printChain();
};


template<class T>
inline double
Metropolis_Hastings_MC<T>::expVal(const double (*f)(T*))
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
Metropolis_Hastings_MC<T>::variance(const double (*f)(T*)){
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
Metropolis_Hastings_MC<T>::printChain(){
	for(unsigned i = 0; i < this->chain.size(); i++)
	{
		std::cout << i << "    " << this->chain[i].logProb() << std::endl;
	}
}


#endif /* Metropolis_Hastings_MC_h */