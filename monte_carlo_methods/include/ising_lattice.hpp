//created by John Stanco 7.2.18


#include "help.hpp"
#include "part_storage.hpp"


#ifndef ISING_LATTICE_H
#define ISING_LATTICE_H


class ising_lattice {
public:
	double E;
	double M;
	double dE;
	double dM;
	std::vector<uint> cluster;
	size_t n_same;
	size_t n_diff;

	virtual void calc_energy() = 0;
	virtual void calc_mag() = 0;
};


class square_ising_lattice : public ising_lattice {
protected:
	arma_tensor<int> 	lat;
	arma_tensor<int>	is_checked;
	void gen_lat();
	void calc_energy();
	void calc_mag();
	void reset() {
		gen_lat();
		calc_energy(); 
		calc_mag(); 
		reset_cluster();
	}
	int is_added(const size_t);
	void add_to_cluster(const size_t);
	void expand_cluster(const size_t, const double);
	void reset_cluster();
public:
	const double h;
	const double J;
	const size_t n_site;
	const size_t dim;
	const arma::uvec r_max;

	square_ising_lattice(const arma::uvec &, double, double);
	void update_obs();
	void flip_spin(const size_t); 
	double flip_and_check(const size_t);
	double nn_energy(const size_t);
	int flip_cluster(const size_t, const double);
	double energy() const { return E; }
	double energy_per_site() const { return E/n_site; }
	double mag_per_site()  const { return M/n_site; } 
	double mag() const { return M; }
	int& operator()(const arma::uvec& indices) { return lat(indices); }
};


#endif /* ISING_LATTICE_H */