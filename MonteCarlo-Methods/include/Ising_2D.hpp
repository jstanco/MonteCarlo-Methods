//created by John Stanco on 5/24/18

#include "world"

#include "world"
#include "help.hpp"
#include "Ising_base.hpp"

#ifndef Ising_2D_h
#define Ising_2D_h

class Ising_2D : public Ising_Base
{
protected:
	int len;
	int wid;
	int nSpins;
	std::vector<int> spins;
	double computeEnergy();
	double computeMag();
public:
	Ising_2D(){}
	Ising_2D(int, double, double, double);
	Ising_2D(int, int, double, double, double);

	virtual int update();
	const std::vector<int>& getSpins() const;
	const int size() const;

	const int& operator()(unsigned, unsigned);
	bool operator==(const Ising_2D &other);
	bool operator >(const Ising_2D&);
};


#endif /* Ising_2D_h */