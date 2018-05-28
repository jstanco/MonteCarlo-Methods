//created by John Stanco on 5/25/18

#include "State_Base_Meat.hpp"

class Ising_Base_Meat : public State_Base_Meat{
	int len;
	double E;
	double M;
	std::vector<int> spins;
	double computeEnergy(double, double);
	double computeMag();
public:
	Ising_Base_Meat(){}
	Ising_Base_Meat(int, double, double);
	Ising_Base_Meat(int, double, double, const Ising_Base_Meat&);
	const std::vector<int>& getSpins() const;

	const int size() const;
	const double energy() const;
	const double mag() const;

	bool operator >(const Ising_Base_Meat&);
	const int& operator[](unsigned);
};