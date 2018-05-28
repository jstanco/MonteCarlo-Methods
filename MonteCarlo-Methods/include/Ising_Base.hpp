//created by John Stanco on 5/23/18

//abstract Ising Model base class

#include "State_Base_Bones.hpp"


#ifndef Ising_Base_h
#define Ising_Base_h


class Ising_Base : public State_Base_Bones
{
protected:
	double T;
	double J;
	double h;
	double E;
	double M;
public:
	Ising_Base(){}
	Ising_Base(double, double, double);

	const double logProb() const;
	virtual const int size() const = 0;
	const double temp() const;
	const double energy() const;
	const double mag() const;

	bool operator ==(const Ising_Base&);
	bool operator >(const Ising_Base&);
	virtual int printData();
};


//--------------------------------------------------------------------------------------//
//									CLASS DEFINITION									//
//--------------------------------------------------------------------------------------//


inline
Ising_Base::Ising_Base(double T, double J, double h) : 
T(T), J(J), h(h), E(0), M(0) {}


inline const double
Ising_Base::energy() const
{
	//printf("%f\n", this->E);
	return this->E;
}


inline const double
Ising_Base::mag() const
{
	return this->M;
}


inline const double
Ising_Base::temp() const
{
	return this->T;
}


inline const double
Ising_Base::logProb() const
{
	return -this->E / (k * this->T);
}


inline int
Ising_Base::printData()
{
	printf("Energy:  %f\nMagnetization:  %f\n", this->energy(), this->mag());
	return 1;
}


#endif /* Ising_Base_h */