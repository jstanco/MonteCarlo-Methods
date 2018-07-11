//created by John Stanco 7.11.18


#include "../include/ising_state.hpp"
#include "../include/metropolis_hastings_mc.hpp"

typedef ising_state<square_ising_lattice> ising_system;
typedef double( * ising_func )( const ising_system & );
typedef std::vector<std::string> str_arr;

template<class T>
double energy( const T & ising )
{
	return ising.energy();
}


template<class T>
double spec_heat( const T & ising )
{
	return ising.spec_heat();
}


template<class T>
double mag( const T & ising )
{
	return ising.mag();
}


template<class T>
double mag_sus( const T & ising )
{
	return ising.mag_sus();
}


std::vector<uint> build_r_max( const size_t dim )
{
	std::vector<uint> 	_r_max( dim );
	for( size_t i = 0; i < dim; i++ )
	{
		_r_max[ i ] = 10;
	}
	return _r_max;
}


template <class T>
int populate_arr( T ( *f )( int ), T * arr, size_t len )
{
	for( size_t i = 0; i < len; i++ )
	{
		arr[ i ] = f( i );
	}
	return 1;
}


int run_T_dependence( 	const size_t iter, 
						const std::vector<uint> & r_max, 
						const str_arr filenames, 
						const ising_func *f_arr, 
						const size_t n_obs )
{

	size_t dim 	= r_max.size();
	double J 	= 1;
	double h 	= 0;

	FILE *pFiles[ n_obs ];
	for( size_t i = 0; i < n_obs; i++ )
	{
		pFiles[ i ] = fopen(filenames[ i ].c_str(), "w.");
	}

	metropolis_hastings_mc<ising_system> MCMC;
	double TMin = .7;
	double TMax = 5;
	size_t Tsteps = 20;
	double mult = (TMax - TMin) / ( Tsteps - 1 );

	for( float i = 0; i < Tsteps; i++ ){
		double T = TMin + i * mult;
		ising_system init( r_max, T, J, h );
		MCMC.run( &init, iter, 5000, 0, false );

		for( size_t j = 0; j < n_obs; j++ )
		{
			fprintf( pFiles[ j ], "%.4f\t\t%.4f\n", T, f_arr[ j ]( init ) );
			//fprintf( pFiles[ j ], "%lf\n", T );
		}	
	}

	for( size_t i = 0; i < n_obs; i++ )
	{
		fclose( pFiles[ i ] );
	}

	return 1;
}



int ising_T_dependence( size_t dim )
{
	
	std::vector<uint> 		 r_max 	= build_r_max(dim);
	size_t 					 iter 	= 1e6;
	str_arr					 filenames( 4 );
	ising_func	 			 f_arr[ 4 ];

	filenames[ 0 ] 	= "../data/Energy_";
	filenames[ 1 ] 	= "../data/Magnetization_";
	filenames[ 2 ] 	= "../data/Specific_Heat_Capacity_";
	filenames[ 3 ] 	= "../data/Magnetic_Susceptibility_";

	f_arr[ 0 ] 		= &energy;
	f_arr[ 1 ] 		= &mag;
	f_arr[ 2 ] 		= &spec_heat;
	f_arr[ 3 ] 		= &mag_sus;

	for( size_t i = 0; i < 4; i++ )
	{
		filenames[ i ] = filenames[ i ] + std::to_string( dim ) + "D.dat";
	}

	run_T_dependence( iter, r_max, filenames, f_arr, 4 );

	return 1;
}


int main(){
	return ising_T_dependence( 2 );
}