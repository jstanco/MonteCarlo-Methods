//created by John Stanco 7.11.18


#include "../../include/help.hpp"
#include "../../../Signal/include/timer.hpp"

int test_2d_loop( size_t nx, size_t ny )
{
	for ( size_t i = 0; i < nx; i++ )
	{
		for ( size_t j = 0; j < ny; j++ )
		{
			size_t k = i * ny + j;
			std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << std::endl;
		}
	}

	std::cout << std::endl;

	for( size_t k = 0; k < nx * ny; k++ )
	{

		size_t j = k % ny;
		size_t i = ( k - j ) / ny;

		std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << std::endl;
	}
	return 1;
}


int test_3d_loop( size_t nx, size_t ny, size_t nz )
{
	for ( size_t i = 0; i < nx; i++ )
	{
		for ( size_t j = 0; j < ny; j++ )
		{
			for ( size_t k = 0; k < nz; k++ )
			{
				size_t l = i * ny * nz + j * nz + k;
				std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tl: " << l << std::endl;
			}
			
		}
	}

	std::cout << std::endl;

	for( size_t l = 0; l < nx * ny * nz; l++ )
	{

		size_t k = l % nz;
		size_t j = ( ( l - k) / nz) % ny;
		size_t i = ( l - j * nz - k ) / ( ny * nz );

		std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tl: " << l << std::endl;
	}
	return 1;
}


template <class T>
T product( const arma::Col<T>& x )
{
	T prod = 1; 
	for( size_t i = 0; i < x.size(); i++ )
	{
		prod *= x( i );
	}
	return prod;
}


arma::uvec compute_stride( const arma::uvec& r_max )
{
  size_t dim = r_max.n_elem;
  arma::uvec stride_len( dim );
  stride_len( dim - 1 ) = 1;
  for( int i = dim - 2; i >= 0; i-- )
  {
    stride_len( i ) = stride_len( i + 1 ) * r_max( i + 1 );
  }
  return stride_len;
}


int to_index( const arma::uvec& indices, const arma::uvec& stride_len )
{
	size_t dim = indices.size();
	size_t ind = 0;
	for ( int i = dim - 1; i >= 0; i-- )
	{
		ind  += indices( i ) * stride_len( i );
	}	
  return ind;
}




int to_index_vec_mult( const arma::uvec& indices, const arma::uvec& stride_len )
{
  return dot( stride_len, indices );
}


arma::uvec to_vec_indices( const arma::uvec& ns, const size_t index )
{
    size_t dim      = ns.n_elem;
    size_t y_i      = 0;
    size_t mult     = 1;
    arma::uvec xs( dim );
    for ( int i = dim - 1; i > 0; i-- )
    {

      xs( i )  = ( ( index - y_i ) / mult ) % ns( i );
      y_i     += xs( i ) * mult;
      mult    *= ns( i );

    }

    xs( 0 ) = ( ( index - y_i ) / mult );
    return xs;
}


//testing indexing over n - dimensional loop
int test_nd_loop( const arma::uvec& ns )
{

	size_t dim 		= ns.size();
	size_t n_tot	= prod( ns );
	size_t loop_max = dim - 1;

	for ( size_t l = 0; l < n_tot; l++ )
	{
    arma::uvec xs = to_vec_indices( ns, l );

		for ( size_t i = 0; i < dim; i++ )
		{
			std::cout << i << ": " << xs( i ) << "\t";
		}

		std::cout << "l: " << l;
		std::cout << "\tl: " << to_index( xs, compute_stride( ns ) ) << std::endl;
	}

	return 1;
}


//random integer indices
arma::uvec rand_indices( const arma::uvec & r_max )
{
  arma::uvec indices( r_max.n_elem );
  for ( size_t i = 0; i < r_max.n_elem; i++ )
  {
    indices( i ) = rand_int( r_max( i ) );
  }
  return indices;
}


int test_index_time( const arma::uvec& r_max )
{
  size_t i1 = 0;
  size_t i2 = 0;
  timer t1, t2;
  arma::uvec stride_len = compute_stride( r_max );
  for( size_t i = 0; i < 100000; i++ )
  {
    arma::uvec indices = rand_indices( r_max );
    t1.start();
    i1 = to_index( indices, stride_len );
    t1.stop();

    t2.start();
    i2 = to_index_vec_mult( indices, stride_len );
    t2.stop();

    if( i1 != i2 ) return 0;
  }

  std::cout << t1.t_sec() / t2.t_sec() << std::endl;

	return 1;
}


int test_loop(){
	size_t d = 2; //d-dimensional lattice
	arma::uvec ns( d );
	for ( size_t i = 0; i < d; i++ )
	{
		ns( i ) = 2;
	}

	return test_nd_loop( ns );
}


int main()
{
	arma::uvec r_max( 3 );
	r_max( 0 ) = 5;
	r_max( 1 ) = 4;
	r_max( 2 ) = 7;

	return test_index_time( r_max );
}