//created by John Stanco on 5/25/18


#include "world"
#include "help.hpp"


#ifndef markov_chain_h
#define markov_chain_h


template<class T>
class markov_chain
{
protected:
	std::vector<T*> chain;
	virtual T* 		transition( T* );
	int 			initialize();
	int 			burn( T* state, uint iter );
	virtual int 	print_data( T* state, uint iter );
public:
	markov_chain(){};
	~markov_chain();
	virtual inline int 	run( T*, uint, uint burn = 0, uint lag = 0, bool save = false, bool print = false );
	virtual int 		clear();

	std::vector<T*>& getChain();
};


template<class T>
inline
markov_chain<T>::~markov_chain(){
	clear();
}


template<class T>
inline T*
markov_chain<T>::transition( T* prev ){
	T* proposed = new T( *prev );
	proposed->update();
	return proposed;
}


template<class T>
inline int
markov_chain<T>::print_data( T* state, uint toGo ){
	printf( "|------------------------------------------------|\n" );
	printf( "Iterations to go:   %d\n", toGo );
	state->print_data();
	return 1;
}


template<class T>
inline int
markov_chain<T>::initialize(){
	clear();
	chain = std::vector<T *>();
	return 1;
}


template<class T>
inline int
markov_chain<T>::burn( T* state, uint _burn ){
	for( uint i = 0; i < _burn; i++ ){
		state->update( true );
	}
	return 1;
}


template<class T>
inline int
markov_chain<T>::run( T *init, uint iter, uint _burn, uint lag, bool save, bool print )
{
	initialize();
	burn( init, _burn );

	for( size_t i = 0; i < iter; i++ )
	{	
		burn( init, lag );
		init->update();
		if( save )
		{
			chain.push_back( new T( *init ) );
		}
		if( print && i % 100000 == 0 )
		{
			print_data( init, iter - i );
		}
	}
	return 1;
}


template<class T>
inline int
markov_chain<T>::clear(){
	for(unsigned i = 0; i < chain.size(); i++)
	{
		if( this->chain[ i ] ) delete this->chain[i];
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