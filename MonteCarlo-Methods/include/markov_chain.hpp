//created by John Stanco on 5/25/18


#include "world"
#include "help.hpp"
#include "Distribution.hpp"


#ifndef markov_chain_h
#define markov_chain_h


template<class T>
class markov_chain
{
protected:
	std::vector<T*> chain;
	virtual T* transition(T*);
	int initialize();
	int burn(T* state, uint iter);
	virtual int printData(T* state, uint iter);
public:
	markov_chain(){};
	~markov_chain();
	virtual inline int run(T*, unsigned, unsigned burn = 0, unsigned lag = 0, bool save = false, bool print = false);
	virtual int clear();

	std::vector<T*>& getChain();
};


template<class T>
inline
markov_chain<T>::~markov_chain(){
	clear();
}


template<class T>
inline T*
markov_chain<T>::transition(T* prev){
	T* proposed = new T(*prev);
	proposed->update();
	return proposed;
}


template<class T>
inline int
markov_chain<T>::printData(T* state, uint toGo){
	printf("|------------------------------------------------|\n");
	printf("Iterations to go:   %d\n", toGo);
	state->printData();
	return 1;
}


template<class T>
inline int
markov_chain<T>::initialize(){
	clear();
	chain = std::vector<T*>();
	return 1;
}


template<class T>
inline int
markov_chain<T>::burn(T* state, uint _burn){
	for(uint i = 0; i < _burn; i++){
		state->update();
	}
	return 1;
}


template<class T>
inline int
markov_chain<T>::run(T *init, unsigned int iter, unsigned int _burn, unsigned int lag, bool save, bool print)
{
	initialize();
	chain.push_back(new T(*init));

	burn(chain[0], _burn);

	if(save){
		for(uint i = 1; i < iter; i++){	
			burn(chain[i - 1], lag);
			chain.push_back(transition(chain[i - 1]));
			if(print && i % 1000000 == 0){
				printData(chain[i], iter - i);
			}
		}
	} else {
		for(uint i = 1; i < iter; i++){	
			burn(chain[0], lag);
			chain[0]->update();
			if(print && i % 1000000 == 0){
				printData(chain[0], iter - i);
			}
		}
	}
	return 1;
}


template<class T>
inline int
markov_chain<T>::clear(){
	for(unsigned i = 0; i < chain.size(); i++)
	{
		if(this->chain[i]) delete this->chain[i];
	}
	this->chain.clear();
	return 1;
}


template<class T>
inline std::vector<T*>&
markov_chain<T>::getChain()
{
	return this->chain;
}
#endif /* markov_chain_h */