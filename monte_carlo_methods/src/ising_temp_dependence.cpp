//created by John Stanco 7.11.18


#include "../include/ising_state.hpp"
#include "../include/markov_chain_mc.hpp"

typedef std::vector<std::string> str_arr;

template<class T>
double energy(const T& ising){
	return mean(ising.energy());
}


template<class T>
double spec_heat(const T& ising){
	return var(ising.spec_heat()) / pow(ising.temp(),2);
}


template<class T>
double mag(const T& ising){
	return ising.mag();
}


template<class T>
double mag_sus(const T & ising){
	return ising.mag_sus();
}


std::vector<uint> build_r_max(const size_t dim){
	std::vector<uint> _r_max(dim);
	for(size_t i = 0; i < dim; i++){
		_r_max[i] = 20;
	}
	return _r_max;
}


str_arr build_filename_string(const std::vector<uint> &rmax){
	str_arr	filenames(4);
	std::string	f_extension = ".dat";
	filenames[0] = "../data/energy";
	filenames[1] = "../data/magnetization";
	filenames[2] = "../data/specific_heat_capacity";
	filenames[3] = "../data/magnetic_susceptibility";
	for(size_t i = 0; i < rmax.size(); i++){
		f_extension = "_" + std::to_string(rmax[i]) + f_extension;
	}
	for(size_t i = 0; i < 4; i++){
		filenames[i] += f_extension;
	}
	return filenames;
}


template<class U>
void run_T_dependence(const size_t iter, const std::vector<uint> &rmax){

	str_arr filenames = build_filename_string(rmax);
	size_t n_obs = 4;
	FILE *pFiles[n_obs];
	double obs[4];

	for(size_t i = 0; i < 4; i++){
		pFiles[i] = fopen(filenames[i].c_str(), "w.");
	}

	markov_chain_mc<U> MCMC;
	double J = 1;
	double h = 0;
	double TMin = .7;
	double TMax = 5;
	size_t Tsteps = 30;
	double mult = (TMax-TMin)/(Tsteps-1);

	for(float i = 0; i < Tsteps; i++){
		double T = TMin + i * mult;
		U init(rmax, T, J, h);
		U rslt = MCMC.run(init, iter, 100, 0, false);

		/// <E> - Energy
		obs[0] = mean(rslt.energy()) / rslt.n_site();
		/// <M> - Magnetization
		obs[1] = mean(rslt.mag()) / rslt.n_site();
		/// Cv - Specific Heat
		obs[2] = var(rslt.energy()) * pow(T,-2) / rslt.n_site();
		/// X - Magnetic Susceptibility
		obs[3] = var(rslt.mag()) * pow(T,-1) / rslt.n_site();

		for(size_t j = 0; j < n_obs; j++){
			fprintf(pFiles[j], "%.4f\t\t%.4f\n", T, obs[j]);
		}	
	}
	for(size_t i = 0; i < n_obs; i++){ fclose(pFiles[i]); }
}


template<class T>
void ising_T_dependence(size_t dim){
	std::vector<uint> r_max = build_r_max(dim);
	size_t 	iter = 3e3;
	run_T_dependence<T>(iter, r_max);
}


int main(){
	timer t;
	t.start();
	ising_T_dependence<square_wolff_ising>(2);
	t.stop();
	t.print();
	return 1;
}