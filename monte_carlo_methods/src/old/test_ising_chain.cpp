//created by John Stanco on 5/12/18

#include "../include/markov_chain_mc.hpp"
#include "../include/ising_1d.hpp"
#include "../include/ising_2d.hpp"
#include <algorithm>

template<class T>
const double
energy(T *chain)
{
	return chain->energy();
}


template<class T>
const double
magnetization(T *chain)
{
	return chain->mag();
}


template<class T>
const double
mag_squared(T *chain)
{
	return pow(chain->mag(), 2);
}


template<class T>
const double
abs_mag(T *chain)
{
	return std::abs(chain->mag());
}


template<class T1>
const double
energyExp(markov_chain_mc<T1>& MCMC, double T)
{
	int size = MCMC.getChain()[0]->size();
	return MCMC.expVal(energy) / size;
}


template<class T1>
const double
magExp(markov_chain_mc<T1>& MCMC, double T)
{
	int size = MCMC.getChain()[0]->size();
	return MCMC.expVal(abs_mag) / size;
}


template<class T1>
const double
specHeatExp(markov_chain_mc<T1>& MCMC, double T)
{
	int size = MCMC.getChain()[0]->size();
	return MCMC.variance(energy) / (k * T * T) / size;
}


template<class T1>
const double
magSusExp(markov_chain_mc<T1>& MCMC, double T)
{
	int size = MCMC.getChain()[0]->size();
	return MCMC.variance(abs_mag) / (k * T) / size;
}


template<class T1>
int
T_Dependence(double L, double J, double h, double TMax, double iter, const double (*f)(markov_chain_mc<T1>&, double), std::string filename)
{
	FILE *pFile = fopen(filename.c_str(), "w.");
	markov_chain_mc<T1> MCMC;
	double TMin = .5;
	double mult = (TMax - TMin) / 100;
	for(float i = 0; i < 100; i++){
		double T = TMin + i * mult;
		T1 *init = new T1(L, T, J, h);
		MCMC.run(init, iter, 0, 0, true, false);
		fprintf(pFile, "%.4f\t\t%.4f\n", T, f(MCMC, T));
		delete init;
	}
	fclose(pFile);
	return 1;
}


template<class T1>
int
runIsingModel(double L, double T, double J, double h, int iter)
{
	markov_chain_mc<T1> MCMC;
	T1 init(L, T, J, h);
	clock_t t = clock();
	MCMC.run(&init, iter, 3e3, 0, true, false);	
	t = clock() - t;		
	printf("\nExpectation value of magnetization of Ising Chain:  %f\n\n", MCMC.expVal(magnetization));
	printf("Expectation value of energy of Ising Chain:  %f\n", MCMC.expVal(energy));
	printf("\nIsing Chain length:  %d\nIterations:  %d\nTime:  %f\n\n", init.size(), iter, ((float)t) / CLOCKS_PER_SEC);
	return 1;
}

//change this such that the updating procedure never copies the state.

//need to compute the energy of a particular state

int main()
{
	//better thing to do would be to pass in an array of function pointers
	//These would all be run on the sate each time it is updated, then
	//This would fall within the MCMC call f(state) on the state each time
	//set seed for random number generator
	srand(time(NULL));

	double L = 10;
	double TMax = 5;
	double J = 1;
	double h = 0;
	double iter = 3e6;

	
	//std::string filename1 = "../data/Energy_1D.dat";
	//std::string filename2 = "../data/Magnetization_1D.dat";
	//std::string filename3 = "../data/Specific_Heat_Capacity_1D.dat";
	//std::string filename4 = "../data/Magnetic_Susceptibility_1D.dat";
	
	//T_Dependence<ising_1d>(L, J, h, TMax, iter, energyExp, filename1);
	//T_Dependence<ising_1d>(L, J, h, TMax, iter, magExp, filename2);
	//T_Dependence<ising_1d>(L, J, h, TMax, iter, specHeatExp, filename3);
	//T_Dependence<ising_1d>(L, J, h, TMax, iter, magSusExp, filename4);
	
	
	
	std::string filename5 = "../data/Energy_2D.dat";
	std::string filename6 = "../data/Magnetization_2D.dat";
	std::string filename7 = "../data/Specific_Heat_Capacity_2D.dat";
	std::string filename8 = "../data/Magnetic_Susceptibility_2D.dat";


	//Design the simulation in such a way that the 
	
	T_Dependence<ising_2d>(L, J, h, TMax, iter, energyExp, filename5);
	T_Dependence<ising_2d>(L, J, h, TMax, iter, magExp, filename6);
	T_Dependence<ising_2d>(L, J, h, TMax, iter, specHeatExp, filename7);
	T_Dependence<ising_2d>(L, J, h, TMax, iter, magSusExp, filename8);
	
	//runIsingModel<ising_2d>(L, 1, J, h, iter);
	return 1;
}