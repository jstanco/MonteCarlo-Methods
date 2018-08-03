//create by John Stanco 7.18.18

#include "../include/ising_lattice.hpp"


square_ising_lattice::square_ising_lattice(const arma::uvec& r_max, double J, double h)
:	r_max(r_max),
	n_site(prod(r_max)),
	dim(r_max.n_elem),
	lat(arma_tensor<int>(r_max)),
	is_checked(arma_tensor<int>(r_max)),
	J(J), h(h)
{ reset(); }


void square_ising_lattice::gen_lat() {
	lat.zeros();
	for (size_t i = 0; i < n_site; i++) {	
		int rand_spin = 2 * rand_int(2) - 1;
		lat(i) = rand_spin;
	}
}


void square_ising_lattice::calc_energy() {
	E = 0; dE = 0;
	for (size_t i = 0; i < n_site; i++) { E += nn_energy(i); }
	E /= 2;
}


void square_ising_lattice::calc_mag() {
	M = 0; dM = 0;
	for (size_t i = 0; i < n_site; i++) { M += lat(i); }
}


void square_ising_lattice::update_obs() {
	E += dE;
	M += dM;
	dE = 0;
	dM = 0;
}


//--------------------- Local-Update ---------------------//


double square_ising_lattice::nn_energy(const size_t index) {
	double E_nn = 0;
	arma::uvec nn_ind	= lat.nn(index);
	for (size_t i = 0; i < 2 * dim; i++) {
		E_nn -= J * lat(nn_ind(i)) * lat(index);
	}
	if (h) { E_nn -= h * lat(index); }
	return E_nn;
}


void square_ising_lattice::flip_spin(const size_t spin_index) {
	lat(spin_index) *= -1;
}


double square_ising_lattice::flip_and_check(const size_t flip_index) {
	flip_spin(flip_index);
	dE = 2 * nn_energy(flip_index);
	if(h) { dE -= 2 * h * lat(flip_index); }
	dM = 2 * lat(flip_index);
	return dE;
}


//--------------- Wolff Cluster Algorithm ----------------//


void square_ising_lattice::reset_cluster() {
	cluster = std::vector<uint>(0);
	n_diff = 0;
	n_same = 0;
	is_checked.zeros();
}


void square_ising_lattice::add_to_cluster(const size_t spin_index) {
	n_same -= is_checked(spin_index);
	is_checked(spin_index) = -1;
	cluster.push_back(spin_index);
}


int square_ising_lattice::is_added(const size_t spin_index) {
	return (is_checked(spin_index) == -1);
}


void square_ising_lattice::expand_cluster(const size_t spin_index, const double flip_prob) {
	add_to_cluster(spin_index);
	arma::uvec nn = lat.nn(spin_index);
	double n_nn = 2 * dim;
	for (size_t i = 0; i < n_nn; i++) {
		if (!is_added(nn(i))) {
			if (lat(nn(i)) == lat(spin_index)) {
				if (flip( flip_prob)){
					expand_cluster(nn(i), flip_prob); 
				} else { n_same++; is_checked(nn(i))++; }
			} else {
				n_diff++;
			}
		}
	}
}


int square_ising_lattice::flip_cluster(const size_t spin_index, const double flip_prob) {
	expand_cluster(spin_index, flip_prob);
	int n_flipped = cluster.size();
	int m_n = n_same - n_diff;

	for (size_t i = 0; i < n_flipped; i++) { flip_spin(cluster[i]); }

	dM = 2 * n_flipped * lat(spin_index);
	dE = 2 * J * m_n;
	dE -= dM * h;
	reset_cluster();
	return m_n;
}

