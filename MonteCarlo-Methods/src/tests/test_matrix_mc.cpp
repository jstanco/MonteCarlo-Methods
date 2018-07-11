//created by John Stanco 7.1.18

#include "../include/matrix_mc_state.hpp"
#include "../include/markov_chain.hpp"


int approx( double x, double y, double etol = 1e-10 )
{
	return ( x == y ) || ( std::abs( x - y ) < etol );
}


int _test_stochastic_matrix()
{
	size_t dim = 3;
	arma::mat st_mat = rand_stochastic_matrix( dim );
	
	for( size_t i = 0; i < dim; i++ )
	{
		double rsum = row_sum( st_mat, i );
		if( !approx( rsum, 1 ) )
		{
			return 0;
		}
	}
	return 1;
}


int test_stochastic_matrix()
{
	for( size_t i = 0; i < 1e5; i++ )
	{
		if( !_test_stochastic_matrix() ){
			return 0;
		}
	}
	return 1;
}


arma::vec norm_vec( arma::vec x )
{
	double sum = vec_sum( x );
	if( sum != 0 ) x /= sum;
	return x;
}


int test_matrix_mc()
{
	size_t dim = 5;
	arma::vec p_0 = arma::randu<arma::vec>( 5 );

	std::cout << vec_sum( p_0 ) << std::endl;
	p_0.print();
	std::cout << "\n" << std::endl;
	p_0 = norm_vec( p_0 );
	std::cout << vec_sum( p_0 ) << std::endl;
	p_0.print();
	std::cout << "\n" << std::endl;
	matrix_mc_state mat_mc( p_0 );

	markov_chain<matrix_mc_state> mc;

	mc.run( &mat_mc, 1e1, 0, 0, false, true );

	return 1;
}


int main()
{
	srand( time( NULL ) );
	//std::cout << test_stochastic_matrix() << std::endl;

	test_matrix_mc();
	

	return 1;
}