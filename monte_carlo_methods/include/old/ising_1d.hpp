//created by John STanco 5/11/18

//implements a 1D Ising chain

#include "world"
#include "help.hpp"
#include "ising_base.hpp"

#ifndef ising_1d_h
#define ising_1d_h

class ising_1d : public ising_base
{
protected:
	int len;
	std::vector<int> spins;
	double computeEnergy();
	double computeMag();
public:
	ising_1d(){}
	ising_1d(int, double, double, double);

	virtual int update();
	const std::vector<int>& getSpins() const;
	const int size() const;

	const int& operator[](unsigned);
	bool operator==(const ising_1d &other);
	bool operator >(const ising_1d&);
};


#endif /* ising_1d_h */


//------------ TODO ------------//
//	Implement 2d & 3D versions
//	Implement potential lattice version