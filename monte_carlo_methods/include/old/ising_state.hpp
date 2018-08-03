// created by John Stanco 7.1.18
// many algorithms defined here: https://arxiv.org/pdf/1404.0209.pdf


#include "state_base.hpp"
#include "ising_lattice.hpp"
#include "markov_chain_mc.hpp"
#include "signal.hpp"


#ifndef ISING_STATE_H
#define ISING_STATE_H


extern double k;


class mcmc_state : public state_base {
protected:
	virtual	bool accept(double);
	virtual double acc_prob(double) = 0;
};

bool
mcmc_state::accept( double acc_prob )
{
	if( acc_prob >= 1 ) { return 1; }
	return flip( acc_prob );
}


template <class U>
class ising_state : public mcmc_state {
	friend class markov_chain_mc<ising_state<U> >;
protected:
	U lat;
	//keeping track of obervables
	Signal::signal M;
	Signal::signal E;
	//arma_tensor G_2; //two-point correlation function
	size_t 	iter;
	size_t 	n_samples;
	double 	T;
	double 	B;
	
	virtual void update_lat();
	virtual void update_obs();
	void init_obs();
	void reset_obs() {}
	void reset_lat() { lat.reset(); update_obs(); }
public:
	ising_state(const std::vector<uint>, const double, const double, const double);
	virtual void set_temp(const double _T) { T = _T; B = pow(1 * T, -1); }
	void prepare(const size_t n_smpl){ M=Signal::signal(n_smpl); E=Signal::signal(n_smpl); iter = 0; n_samples=n_smpl; }
	/// burn discards move 
	void update(const bool burn = false); 
	void reset();
	/// observables
	double n_site() const { return lat.n_site; }
	Signal::signal mag() const { return M; }
	Signal::signal energy() const { return E; }
	void print_data() const;
};


/// Metropolis-Hastings
template<class U>
class ising_state_mh : public ising_state<U> {
protected:
	double acc_prob(double dE) {
		return exp(-this->B * dE);
	}
public:
	ising_state_mh(const std::vector<uint> r_max, 
						const double T, 
						const double J, 
						const double h) : 
	ising_state<U>(r_max, T, J, h) {}
};


template<class U>
class ising_state_glauber : public ising_state<U> {
protected:
	double acc_prob(double dE) {
		double b = exp(-this->B * dE);
		return b / (1 + b);
	}
public:
	ising_state_glauber(const std::vector<uint> r_max, 
						const double T, 
						const double J, 
						const double h) : 
	ising_state<U>(r_max, T, J, h) {}
};


template<class U>
class ising_state_heat_bath : public ising_state<U> {
protected:
	double acc_prob(double E_nn) {
		double dS = -this->B * 2 * E_nn;
		double b = exp(this->B * dS);
		return b / (1 + b);
	}
public:
	ising_state_heat_bath(const std::vector<uint> r_max, 
						const double T, 
						const double J, 
						const double h) : 
	ising_state<U>(r_max, T, J, h) {}

	void update_lat() {
		int flip_index = rand_int(this->lat.n_site);
		double E_nn = this->lat.nn_energy(flip_index);
		if (!this->accept(acc_prob(E_nn))) {
			this->lat.flip_and_check(flip_index);
			this->lat.update_obs();
		}
	}
};


//wolff cluster algorithm
template<class U>
class ising_state_wolff : public ising_state_mh<U> {													
protected:
	double flip_prob;
	double calc_flip_prob(const double _T, const double J){
		return 1 - exp(-2 * J / (1 * _T));
	}
public:												
	ising_state_wolff(const std::vector<uint> r_max, 
						const double T, 
						const double J, 
						const double h ) : 
	flip_prob( calc_flip_prob(T, J)),
	ising_state_mh<U>(r_max, T, J, h) {}

	void set_temp(const double _T){
		this->T = _T;
		this->B = pow(1 * _T, -1);
		flip_prob = calc_flip_prob(this->lat.J, _T);
	}

	void update_lat() {
		int flip_index = rand_int(this->lat.n_site);
		this->lat.flip_cluster(flip_index, flip_prob);
		this->lat.update_obs();
	}
};


template <class U>
ising_state<U>::ising_state(const std::vector<uint> r_max, 
							const double T, 
							const double J, 
							const double h)
:	lat(U(to_arma(r_max), J, h)),
	T(T),
	B(pow(1 * T, -1)),
	iter(0) {}


template<class T>
void ising_state<T>::update_obs() {
	M.record(std::abs(lat.mag()), iter);
	//M.record(lat.mag_per_site(), iter);
	E.record(lat.energy(), iter);
}


template <class T>
void ising_state<T>::update_lat() {
	double flip_index = rand_int(lat.n_site);
	double dE = lat.flip_and_check(flip_index);
	if (accept(acc_prob(dE))) { lat.update_obs(); } 
	else { lat.flip_spin(flip_index); }
}


template <class T>
void ising_state<T>::update(const bool burn) {
	update_lat();
	if (!burn) {
		update_obs();
		iter++;
	}
}


template <class T>
void ising_state<T>::reset() {
	reset_obs();
	reset_lat();
	M = Signal::signal();
	E = Signal::signal();
	iter = 0;
}


template<class T>
void ising_state<T>::print_data() const {
	printf("Energy:  %lf\nMagnetization:  %lf\n", lat.energy(), lat.mag());
}


typedef ising_state_mh<square_ising_lattice> 		 square_mh_ising;
typedef ising_state_glauber<square_ising_lattice> 	 square_glauber_ising;
typedef ising_state_heat_bath<square_ising_lattice>  square_heat_bath_ising;
typedef ising_state_wolff<square_ising_lattice>		 square_wolff_ising;

#endif /* ISING_STATE_H */