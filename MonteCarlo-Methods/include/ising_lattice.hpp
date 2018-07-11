//created by John Stanco 7.2.18


#include "help.hpp"


#ifndef ising_lattice_hpp
#define ising_lattice_hpp


template<class T>
arma::Col<T> to_arma( const std::vector<T> & x )
{
	arma::Col<T> y( x.size() );
	for( size_t i = 0; i < x.size(); i++ )
	{
		y( i ) = x[ i ];
	}
	return y;
}


template<class T>
class arma_tensor
{
private:
	size_t 	dim;
	size_t	n_elem;
	arma::Col<T> 	data;
	arma::uvec		n_max;
	arma::uvec		stride_len;

	int compute_stride()
	{
		stride_len( dim - 1 ) = 1;
		for( int i = dim - 2; i >= 0; i-- )
		{
			stride_len( i ) = stride_len( i + 1 ) * n_max( i + 1 );
		}
		return 1;
	}

	size_t to_index( const arma::uvec & indices )
	{
		return dot( indices, stride_len );
	}

	arma::uvec to_vec_indices( const size_t index )
	{
	    size_t 		y_i = 0;
	    arma::uvec 	xs( dim );

	    for ( int i = dim - 1; i > 0; i-- )
	    {
	     	xs( i )  = ( ( index - y_i ) / stride_len( i ) ) % n_max( i );
	     	y_i     += xs( i ) * stride_len( i );
	    }

	    xs( 0 ) = ( ( index - y_i ) / stride_len( 0 ) );
	    return xs;
	}	


public:

	arma_tensor( const arma::uvec & n_max ) 
	:	dim( n_max.n_elem 			  		) ,
		n_elem( prod( n_max ) 		 		) ,
		data( arma::Col<T>( prod( n_max ) ) ) ,
		n_max( n_max 				  		) ,
		stride_len( arma::uvec( n_max.n_elem ) )
	{ compute_stride(); }

	T& operator()( const arma::uvec & indices )
	{
		return data( to_index ( indices ) );
	}

	T& operator()( const size_t index ){
		return data( index );
	}

	arma::uvec nn( const size_t index ) //returns nearest neighbors of cell at particular index;
	{
		arma::uvec _nn( 2 * dim );
		arma::uvec indices = to_vec_indices( index );

		for( size_t i = 0; i < 2 * dim; i += 2 )
		{
			size_t j 	 = i / 2;
			size_t k 	 = indices( j );
			_nn( i ) 	 = index + ( ( k + 1 ) % n_max ( j ) - k ) 				* stride_len( j );
			_nn( i + 1 ) = index + ( ( k + n_max( j ) - 1 ) % n_max( j ) - k ) 	* stride_len( j );
		}
		return _nn;
	}

};


class ising_lattice
{
public:

	double E;
	double M;
	double dE;
	double dM;

	virtual int calc_energy() = 0;
	virtual int calc_mag() = 0;

	virtual double flip_rand() = 0; //returns change in action
	virtual int rev_flip() = 0;	//reverses last flip
};


//-------------------------------------------------//



//-------------------------------------------------//


class square_ising_lattice : public ising_lattice
{
private:
	const size_t 		n_site;
	const size_t 		dim;
	const arma::uvec 	r_max;
	arma_tensor<int> 	lat;
	double 				h;
	double 				J;
	size_t 				flip_index; //saves last move in case of reversal.
	
	int 	gen_lat();
	double 	nn_energy( const size_t );
	int		calc_energy();
	int 	calc_mag();

	int 	reset();
	
public:

	square_ising_lattice( const arma::uvec &, double, double );

	int 	update_obs();
	double 	flip_rand();
	int    	rev_flip(); 

	int&	operator()( const arma::uvec & );
	double 	energy() 	const;
	double 	mag() 		const;
};


square_ising_lattice::square_ising_lattice( const arma::uvec & r_max , double J, double h )
:	r_max( r_max 			 		) ,
	n_site( prod( r_max )    		) ,
	dim( r_max.n_elem		 		) ,
	lat( arma_tensor<int> ( r_max ) ) ,
	J( J 							) ,
	h( h 							)
{ 
	gen_lat(); 
	calc_energy(); 
	calc_mag(); 
}


int
square_ising_lattice::gen_lat()
{
	for( size_t i = 0; i < n_site; i++ )
	{	
		int rand_spin = 2 * rand_int( 2 ) - 1;
		lat( i ) 	  = rand_spin;
	}
	return 1;
}


int
square_ising_lattice::calc_energy()
{
	E = 0;
	for( size_t i = 0; i < n_site; i++ )
	{
		E += nn_energy( i );
	}
	E /= 2;
	return 1;
}


int
square_ising_lattice::calc_mag()
{
	M = 0;
	for( size_t i = 0; i < n_site; i++ )
	{
		M += lat( i );
	}
	return 1;
}


int
square_ising_lattice::update_obs()
{
	
	E += dE;
	M += dM;
	//can have whole list of obs.
	return 1;
}


double
square_ising_lattice::nn_energy( const size_t index )
{
	double E_nn 		= 0;
	arma::uvec nn_ind 	= lat.nn( index );
	for( size_t i = 0; i < 2 * dim; i++ )
	{
		E_nn -= J * lat( nn_ind( i ) ) * lat( index );
	}
	return E_nn;
}


double
square_ising_lattice::flip_rand()
{
	flip_index 			= rand_int( n_site );
	lat( flip_index )  *= -1;

	dE  = 2 * nn_energy( flip_index );
	dE -= 2 * h * lat( flip_index );
	dM  = 2 * lat( flip_index );

	return dE;

}


int
square_ising_lattice::rev_flip()
{
	lat( flip_index ) *= -1;
	return 0;
}


double
square_ising_lattice::energy() const
{
	return E;
}


double
square_ising_lattice::mag() const
{
	return M;
}


int&
square_ising_lattice::operator()( const arma::uvec & indices )
{
    return lat( indices );
}


int
square_ising_lattice::reset()
{
	flip_index = 0;
	E 	= 0;
	M 	= 0;
	dE 	= 0;
	dM 	= 0;
	lat = arma_tensor<int>( r_max );
	return 1;
}


#endif /* ising_lattice_hpp */