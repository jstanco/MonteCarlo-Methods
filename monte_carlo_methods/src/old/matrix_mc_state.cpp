//created by John Stanco 7.1.18

#include "../include/matrix_mc_state.hpp"

double
vec_sum( const arma::vec& x)
{
	double sum = 0;
	for( size_t i = 0; i < x.n_elem; i++ )
	{
		sum += x( i );
	}
	return sum;
}


double 
row_sum( const arma::mat& A, const size_t row )
{
	double rsum  = 0;
	for( size_t col = 0; col < A.n_cols; col++ )
	{
		rsum += A( row, col );
	}
	return rsum;
}


arma::mat
rand_stochastic_matrix( size_t dim )
{
	arma::mat st_mat = arma::randu<arma::mat>( dim, dim );
	for( size_t i = 0; i < dim; i++ )
	{
		double rsum = row_sum( st_mat, i );
		st_mat.row( i ) /= rsum;
	}
	return st_mat;
}


matrix_mc_state::matrix_mc_state( const arma::vec& p_0 ) : 
	tn_mat( rand_stochastic_matrix( p_0.n_elem ) ) ,
	probs( p_0 ) {}


matrix_mc_state::matrix_mc_state( const arma::mat& tn_mat, const arma::vec& p_0 ) : 
	tn_mat( check_stochastic( tn_mat ) ) ,
	probs( p_0 ) {}


arma::mat
matrix_mc_state::check_stochastic( const arma::mat& A )
{
	//first check for row sum...
	if( A.n_rows != A.n_cols )
	{
		throw "matrix_mc_state::check_stochastic : stochastic matrix must be square!";
	}
	for( size_t row = 0; row < A.n_rows; row++ )
	{
		double rsum = 0;
		for( size_t col = 0; col < A.n_cols; col++ )
		{
			if( A( row, col ) < 0 )
			{
				throw "matrix_mc_state::check_stochastic : stochastic matrix must be positive!";
			}
		 	rsum += A( row, col );
		}

		if( rsum != 1 )
		{
			throw "matrix_mc_state::check_stochastic : stochastic matrix rows must sum to one!";
		}
	}

	return A;
}


int
matrix_mc_state::update()
{
	probs = tn_mat * probs;
	std::cout << vec_sum( probs ) << std::endl;
	probs.print();
	std::cout << std::endl;
	return 1;
}


int
matrix_mc_state::print_data()
{
	std::cout << "\nProbabilities:\n" << std::endl;
	probs.print();
	return 1;
}


int
matrix_mc_state::print_mat()
{
	std::cout << "Transition Matrix:" << std::endl;
	tn_mat.print();
	return 1;
}