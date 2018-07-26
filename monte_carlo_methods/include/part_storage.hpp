//created by John Stanco 7.13.18
// Uses armadillo vectors to represent tensor (rank n)

#ifndef ARMA_TENSOR
#define ARMA_TENSOR

#include "mc_methods.hpp"

template<class T>
arma::Col<T> to_arma(const std::vector<T>& x){
	arma::Col<T> y(x.size());
	for(size_t i = 0; i < x.size(); i++)
	{
		y(i) = x[i];
	}
	return y;
}


template<class T>
class arma_tensor{
private:
	size_t dim;
	size_t n_elem;
	arma::Col<T> data;
	arma::uvec n_max;
	arma::uvec stride_len;

	void compute_stride(){
		stride_len(dim - 1) = 1;
		for(int i = dim - 2; i >= 0; i--){
			stride_len(i) = stride_len(i + 1) * n_max(i + 1);
		}
	}

	size_t to_index(const arma::uvec & indices){ return dot(indices, stride_len); }
	arma::uvec to_vec_indices(const size_t index){
	    size_t y_i = 0;
	    arma::uvec xs(dim);
	    for(int i = dim - 1; i > 0; i--){
	     	xs(i) = ((index - y_i) / stride_len(i)) % n_max(i);
	     	y_i += xs(i) * stride_len(i);
	    }
	    xs(0) = ((index - y_i) / stride_len(0));
	    return xs;
	}	
public:
	arma_tensor(const arma::uvec & n_max) 
	:	dim(n_max.n_elem) ,
		n_elem(prod(n_max)),
		data(arma::Col<T>(prod(n_max))),
		n_max(n_max),
		stride_len(arma::uvec(n_max.n_elem))
	{ compute_stride(); }
	T& operator()(const arma::uvec& indices){ return data(to_index(indices)); }
	T& operator()(const size_t index){ return data(index); }
	/// returns nearest neighbors of cell at particular index;
	arma::uvec nn(const size_t index){
		arma::uvec _nn(2 * dim);
		arma::uvec indices = to_vec_indices(index);
		for(size_t i = 0; i < 2 * dim; i += 2){
			size_t j = i / 2;
			size_t k = indices(j);
			_nn(i) = index + ((k + 1) % n_max(j) - k) * stride_len(j);
			_nn(i + 1) = index + ((k + n_max(j) - 1) % n_max(j) - k)*stride_len(j);
		}
		return _nn;
	}
	void zeros(){ data.zeros(); }
};

#endif /* ARMA_TENSOR */