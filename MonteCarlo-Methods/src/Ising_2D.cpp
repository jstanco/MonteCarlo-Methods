//created by John Stanco on 5/26/18

#include "../include/ising_2d.hpp"

#define SUB operator()


//------------------------------------------------------------------------//
//	  							   Helpers								  //
//------------------------------------------------------------------------//


std::vector<int>
random_ising_2d(unsigned int len, unsigned int wid)
{
	unsigned int size = len * wid;
	std::vector<int> chain(size);
	for(unsigned i = 0; i < size; i++)
	{
		int randInt = 2 * (rand() % 2) - 1;
		chain[i] = randInt;
	}
	return chain;
}


std::vector<int>
ordered_ising_2d(unsigned int len, unsigned int wid)
{
	unsigned int size = len * wid;
	std::vector<int> chain(size);
	for(unsigned i = 0; i < size; i++)
	{
		chain[i] = 1;
	}
	return chain;
}


//------------------------------------------------------------------------//
//	  							   Methods								  //
//------------------------------------------------------------------------//


//square <--> len = wid
ising_2d::ising_2d(int size, double T, double J, double h) : 
len(size), wid(size), nSpins(len * wid), ising_base(T, J, h)
{
	this->spins = random_ising_2d(len, wid);
	this->E = this->computeEnergy();
	this->M = this->computeMag();
}


ising_2d::ising_2d(int len, int wid, double T, double J, double h) : 
len(len), wid(wid), nSpins(len * wid), ising_base(T, J, h)
{
	this->spins = random_ising_2d(len, wid);
	this->E = this->computeEnergy();
	this->M = this->computeMag();
}


int
ising_2d::update(){
	unsigned int i = rand() % len;
	unsigned int j = rand() % wid;
	this->spins[i * len + j] *= -1;
	double dM = 2 * this->SUB(i, j);
	double dE = 0;
	if(this->len < 2){
		dE = 0;
	} else {
		dE -= 2 * J * this->SUB(i, j) * this->SUB((i + this->len - 1) % this->len, j);
		dE -= 2 * J * this->SUB(i, j) * this->SUB((i + 1) % this->len, j);
		dE -= 2 * J * this->SUB(i, j) * this->SUB(i, (j + this->wid - 1) % this->wid);
		dE -= 2 * J * this->SUB(i, j) * this->SUB(i, (j + 1) % this->wid);
		dE -= 2 * h * this->SUB(i, j);
	}
	if(!this->accept(exp(-dE * this->B))){
		this->spins[i * len + j] *= -1;
	} else {
		this->E += dE;
		this->M += dM;
	}
	return 1;
}


const int
ising_2d::size() const
{
	return this->nSpins;
}


inline const std::vector<int>&
ising_2d::getSpins() const
{
	return this->spins;
}


double
ising_2d::computeEnergy()
{
	double Hs = 0;
	double Hm = 0;
	for(unsigned i = 0; i < this->len; i++){
		for(unsigned j = 0; j < this->wid; j++){
			Hs += this->SUB(i, j) * this->SUB((i + 1) % this->len, j);
			Hs += this->SUB(i, j) * this->SUB(i, (j + 1) % this->wid);
			Hm += this->SUB(i, j);
		}
	}
	return E = -Hs * this->J - Hm * this->h;
}


double
ising_2d::computeMag(){
	double acc = 0;
	for(unsigned i = 0; i < this->nSpins; i++)
	{
		acc += this->spins[i];
	}
	return acc;
}


const int&
ising_2d::operator()(unsigned int i, unsigned int j)
{
	if((i > this->len - 1) || (j > this->wid - 1))
	{
		throw "|  function: ising_2d::operator ()  |  file: ising_2d.cpp  |  error: Index out of bounds |";
	}
	return this->spins[i * this->len + j];
}


inline bool
ising_2d::operator==(const ising_2d &other)
{
	if((this->len != other.len) || (this->wid != other.wid))
	{
		return false;
	}
	return this->spins == other.spins;
}