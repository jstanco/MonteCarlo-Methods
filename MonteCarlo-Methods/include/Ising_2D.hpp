//created by John Stanco on 5/24/18

#include "world"

#include "world"
#include "help.hpp"
#include "ising_base.hpp"

#ifndef ising_2d_h
#define ising_2d_h

class ising_2d : public ising_base
{
protected:
	int len;
	int wid;
	int nSpins;
	std::vector<int> spins;
	double computeEnergy();
	double computeMag();
public:
	ising_2d(){}
	ising_2d(int, double, double, double);
	ising_2d(int, int, double, double, double);

	virtual int update();
	const std::vector<int>& getSpins() const;
	const int size() const;

	const int& operator()(unsigned, unsigned);
	bool operator==(const ising_2d &other);
	bool operator >(const ising_2d&);
};


#endif /* ising_2d_h */