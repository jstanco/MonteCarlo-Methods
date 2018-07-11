//created by John Stanco 7.1.18

#include "state_base.hpp"
#include "ising_lattice.hpp"


#ifndef ising_state_hpp
#define ising_state_hpp


extern double k;


class mcmc_state : public state_base
{
protected:
	bool 	accept( double );
public:

};

bool
mcmc_state::accept( double logProb )
{
	if( logProb >= 0 ) { return 1; }
	return flip( exp( logProb ) );
}


//-------------------------------------------------//


template <class T1>
class ising_state : public mcmc_state
{
private:

	T1 lat;

	//keeping track of obervables
	arma::mat obs;
	size_t 	iter;

	double 	T;
	double 	B;
	
	int 	update_lat();
	int 	update_moments();
	int 	update_obs();

	int 	init_obs();
	int 	reset_obs();
	int 	reset_lat();

public:
			ising_state( const std::vector<uint>, 
						const double J, 
						const double h, 
						const double T );

	int 	set_temp( const double );
	int 	update( bool burn = false ); //burn discards move  

	int 	reset();

	//observables
	double 	energy() 	const;
	double 	spec_heat() const;
	double	mag_sus() 	const;
	double 	mag() 		const;

	int 	print_data();

};


template <class T1>
ising_state<T1>::ising_state( const std::vector<uint> r_max, 
							const double T, 
							const double J, 
							const double h )
:	lat 	( T1( to_arma( r_max ), J, h ) ) ,
	obs 	( arma::mat( 2, 3 ) ) ,
	T 		( T 				) ,
	B 		( pow( 1 * T, -1 ) 	) ,
	iter	( 0					) 
{ obs.zeros(); }	//calculates lattice energy -> sets observables;


template<class T>
int
ising_state<T>::update_obs()
{
	obs( 0, 0 ) = lat.energy();
	obs( 1, 0 ) = std::abs( lat.mag() );
	return 1;
}


template <class T>
int
ising_state<T>::update_lat()
{
	double dS = B * lat.flip_rand();
	if( accept( -dS ) ){
		lat.update_obs();
		return 1;
	}
	return lat.rev_flip();
}


template <class T>
int
ising_state<T>::update_moments()
{
	for( size_t i = 0; i < obs.n_rows; i++ )
	{
		//updating moments
		for( size_t j = 1; j < obs.n_cols; j++ )
		{
			obs( i, j ) = ( obs( i, j ) * iter + pow( obs( i, 0 ), j ) ) / (float)( iter + 1 );
		}
	}
	return 1;
}


template <class T>
int
ising_state<T>::update( bool burn )
{
	update_lat();
	update_obs();
	if( !burn )
	{
		update_moments();
		iter++;
	}
	return 1;
}


template <class T>
int
ising_state<T>::reset_obs()
{
	obs.zeros();
}


template <class T>
int
ising_state<T>::reset_lat()
{
	lat.reset();
	update_obs();
}


template <class T>
int
ising_state<T>::reset()
{
	reset_obs();
	reset_lat();
	iter = 0;
	return 1;
}


template <class T1>
int
ising_state<T1>::set_temp( const double _T )
{
	T = _T;
	B = pow( 1 * T, -1 );
	return 1;
}


template <class T>
double
ising_state<T>::energy() const
{
	return obs( 0, 1 );
}


template <class T>
double
ising_state<T>::spec_heat() const
{
	//can store this value
	return ( obs( 0, 2 ) - pow( energy(), 2 ) ) / ( T * T );
}


template <class T>
double
ising_state<T>::mag() const
{
	return obs( 1, 1 );
}


template <class T>
double
ising_state<T>::mag_sus() const
{
	//can store this value
	return ( obs( 1, 2 ) - pow( mag(), 2 ) ) / T; 
}


template<class T>
int
ising_state<T>::print_data()
{
	printf("Energy:  %f\nMagnetization:  %f\n", energy(), mag());
	return 1;
}


#endif /* ising_state_hpp */