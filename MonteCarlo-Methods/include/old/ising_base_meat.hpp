//created by John Stanco on 5/25/18

#include "State_Base_Meat.hpp"

class ising_base_Meat : public State_Base_Meat{
	int len;
	double E;
	double M;
	std::vector<int> spins;
	double computeEnergy(double, double);
	double computeMag();
public:
	ising_base_Meat(){}
	ising_base_Meat(int, double, double);
	ising_base_Meat(int, double, double, const ising_base_Meat&);
	const std::vector<int>& getSpins() const;

	const int size() const;
	const double energy() const;
	const double mag() const;

	bool operator >(const ising_base_Meat&);
	const int& operator[](unsigned);
};