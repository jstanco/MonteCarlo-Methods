//created by John Stanco on 5/25/18


#include "world"
#include "help.hpp"
#include "Distribution.hpp"


#ifndef MarkovChain_h
#define MarkovChain_h


template<class T>
class MarkovChain
{
protected:
	std::vector<T*> chain;
	virtual T* transition(T*);
public:
	MarkovChain(){};
	~MarkovChain();
	virtual inline int run(T*, unsigned, unsigned burn = 0, unsigned lag = 0, bool print = false);
	virtual int clear();

	std::vector<T*>& getChain();
};


template<class T>
inline
MarkovChain<T>::~MarkovChain()
{
	this->clear();
}


template<class T>
inline T*
MarkovChain<T>::transition(T* prev)
{
	return new T(*prev, true);
}


template<class T>
inline int
MarkovChain<T>::run(T *init, unsigned int iterations, unsigned int burn, unsigned int lag, bool print)
{
	this->clear();
	this->chain = std::vector<T*>();
	this->chain.push_back(new T(*init));
	for(unsigned i = 1; i < iterations; i++)
	{
		this->chain.push_back(this->transition(this->chain[i - 1]));
		if(print && i % 1000000 == 0)
		{
			printf("|------------------------------------------------|\n");
			printf("Iterations to go:   %d\n", iterations - i);
			this->chain[i]->printData();
		}
	}
	return 1;
}


template<class T>
inline int
MarkovChain<T>::clear(){
	for(unsigned i = 0; i < chain.size(); i++)
	{
		if(this->chain[i]) delete this->chain[i];
	}
	this->chain.clear();
	return 1;
}


template<class T>
inline std::vector<T*>&
MarkovChain<T>::getChain()
{
	return this->chain;
}
#endif /* MarkovChain_h */