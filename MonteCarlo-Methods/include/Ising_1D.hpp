//created by John STanco 5/11/18

//implements a 1D Ising chain

#include "world"
#include "help.hpp"
#include "Ising_base.hpp"

#ifndef Ising_1D_h
#define Ising_1D_h

class Ising_1D : public Ising_Base
{
protected:
	int len;
	std::vector<int> spins;
	double computeEnergy();
	double computeMag();
public:
	Ising_1D(){}
	Ising_1D(int, double, double, double);
	Ising_1D(const Ising_1D&, bool);

	const std::vector<int>& getSpins() const;
	const int size() const;

	const int& operator[](unsigned);
	bool operator==(const Ising_1D &other);
	bool operator >(const Ising_1D&);
};


#endif /* Ising_1D_h */


//------------ TODO ------------//
//	Implement 2d & 3D versions
//	Implement potential lattice version