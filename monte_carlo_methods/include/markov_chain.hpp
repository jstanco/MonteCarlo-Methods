//created by John Stanco on 5/25/18


#include "mc_methods.hpp"
#include "help.hpp"

#ifndef MARKOV_CHAIN_h
#define MARKOV_CHAIN_h


template<class T>
class markov_chain
{
protected:
	std::vector<T*> chain;
	virtual T* transition(const T&);
	void initialize();
	void burn(T& state, size_t iter);
	virtual void print_data(const T& state, size_t iter);
public:
	markov_chain(){};
	~markov_chain();
	virtual T run(T init, size_t iter, size_t burn = 0, size_t lag = 0, bool save = false, bool print = false );
	virtual void clear();

	std::vector<T*>& getChain();
};


template<class T>
inline
markov_chain<T>::~markov_chain(){
	clear();
}


template<class T>
inline T*
markov_chain<T>::transition(const T& prev){
	T* proposed = new T(prev);
	proposed->update();
	return proposed;
}


template<class T>
inline void
markov_chain<T>::print_data(const T& state, size_t toGo){
	printf("|------------------------------------------------|\n");
	printf("Iterations to go:   %lu\n", toGo);
	state.print_data();
}


template<class T>
inline void
markov_chain<T>::initialize(){
	clear();
	chain = std::vector<T *>();
}


template<class T>
inline void
markov_chain<T>::burn(T& state, size_t n){
	for(size_t i = 0; i < n; i++){
		state.update(true);
	}
}


template<class T>
inline T
markov_chain<T>::run(T init, size_t iter, size_t _burn, size_t lag, bool save, bool print){
	initialize();
	init.prepare(iter);
	burn(init, _burn);

	for(size_t i = 0; i < iter; i++){	
		burn(init, lag);
		init.update();
		if(save){
			chain.push_back(new T(init));
		}
		if(print && i % 100000 == 0){
			print_data(init, iter - i);
		}
	}
	return init;
}


template<class T>
inline void
markov_chain<T>::clear(){
	for(size_t i = 0; i < chain.size(); i++){
		if(this->chain[i]) delete this->chain[i];
	}
	this->chain.clear();
}


template<class T>
inline std::vector<T*>&
markov_chain<T>::getChain(){
	return this->chain;
}
#endif /* MARKOV_CHAIN_h */