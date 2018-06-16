//created by John Stanco on 5/24/18

#include "../include/ising_1d.hpp"

double k = 1;


//------------------------------------------------------------------------//
//	  							   Helpers								  //
//------------------------------------------------------------------------//


std::vector<int>
random_ising_1d(unsigned int len)
{
	std::vector<int> chain(len);
	for(unsigned i = 0; i < len; i++){
		int randInt = 2 * (rand() % 2) - 1;
		chain[i] = randInt;
	}
	return chain;
}


std::vector<int>
ordered_ising_1d(unsigned int len)
{
	std::vector<int> chain(len);
	for(unsigned i = 0; i < len; i++){
		chain[i] = 1;
	}
	return chain;
}


//------------------------------------------------------------------------//
//	  							   Methods								  //
//------------------------------------------------------------------------//


ising_1d::ising_1d(int len, double T, double J, double h) : 
len(len), ising_base(T, J, h)
{
	this->spins = random_ising_1d(len);
	this->E = this->computeEnergy();
	this->M = this->computeMag();
}


int
ising_1d::update(){
	unsigned int index = rand() % this->len;
	//try change
	this->spins[index] *= -1;
	double dM = 2 * this->spins[index];
	double dE = 0;

	if(this->len < 2){
		dE = 0;
	} else {
		dE -= 2 * J * this->spins[index] * this->spins[(index + this->len - 1) % this->len];
		dE -= 2 * J * this->spins[index] * this->spins[(index + 1) % this->len];
		dE -= 2 * h * this->spins[index];
	}
	if(!this->accept(exp(-dE * this->B))){
		spins[index] *= -1;
	} else {
		this->E += dE;
		this->M += dM;
	}
	return 1;
}


const int
ising_1d::size() const
{
	return this->len;
}


inline const std::vector<int>&
ising_1d::getSpins() const
{
	return this->spins;
}


double
ising_1d::computeEnergy()
{
	double Hs = spins[this->len - 1] * spins[0];
	double Hm = this->spins[0];
	for(unsigned i = 1; i < this->len; i++)
	{
		Hs += this->spins[i - 1] * this->spins[i];
		Hm += this->spins[i];
	}
	return E = -Hs * this->J - Hm * this->h;
}


double
ising_1d::computeMag(){
	double acc = 0;
	for(unsigned i = 0; i < this->len; i++)
	{
		acc += this->spins[i];
	}
	return acc;
}


const int&
ising_1d::operator[](unsigned int i)
{
	if(i > this->len - 1)
	{
		throw "|  function: ising_1d::operator []  |  file: ising_1d.cpp  |  error: Index out of bounds |";
	}
	return this->spins[i];
}


inline bool
ising_1d::operator==(const ising_1d &other)
{
	if(this->len != other.len)
	{
		return false;
	}
	return this->spins == other.spins;
}