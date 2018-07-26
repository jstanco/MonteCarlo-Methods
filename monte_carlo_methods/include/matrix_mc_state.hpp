//created by John Stanco 7.1.18

#include "help.hpp"
#include "state_base.hpp"

#ifndef matrix_mc_state_hpp
#define matrix_mc_state_hpp


double		vec_sum(const arma::vec& x);
double 		row_sum(const arma::mat& A, const size_t row );
arma::mat rand_stochastic_matrix(size_t dim);



class matrix_mc_state : public state_base
{
private:
	const arma::mat tn_mat;
	arma::vec 		probs;

	arma::mat check_stochastic(const arma::mat&);

public:
	matrix_mc_state(const arma::vec&);
	matrix_mc_state(const arma::mat&, const arma::vec&);

	//int burn( size_t );
	int update();
	int print_data();
	int print_mat();
};


#endif /* matrix_mc_state_hpp */
