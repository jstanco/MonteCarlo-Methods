//created by John Stanco on 5/23/18

//abstract Ising Model base class

#include "state_base.hpp"


#ifndef ising_base_h
#define ising_base_h


class ising_base : public state_base
{
protected:
	double T;
	double B;
	double J;
	double h;
	double E;
	double M;
public:
	ising_base(){}
	ising_base(double, double, double);

	virtual const int size() const = 0;
	const double temp() const;
	const double energy() const;
	const double mag() const;

	bool operator ==(const ising_base&);
	bool operator >(const ising_base&);
	virtual int printData();
};


//--------------------------------------------------------------------------------------//
//									CLASS DEFINITION									//
//--------------------------------------------------------------------------------------//


inline
ising_base::ising_base(double T, double J, double h) : 
T(T), J(J), h(h), E(0), M(0), B(pow(k * T, -1)) {}


inline const double
ising_base::energy() const
{
	//printf("%f\n", this->E);
	return this->E;
}


inline const double
ising_base::mag() const
{
	return this->M;
}


inline const double
ising_base::temp() const
{
	return this->T;
}


inline int
ising_base::printData()
{
	printf("Energy:  %f\nMagnetization:  %f\n", this->energy(), this->mag());
	return 1;
}


#endif /* ising_base_h */