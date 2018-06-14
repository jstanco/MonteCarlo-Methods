//created by John Stanco on 5/26/18

#include "../include/Ising_2D.hpp"

#define SUB operator()


//------------------------------------------------------------------------//
//	  							   Helpers								  //
//------------------------------------------------------------------------//


std::vector<int>
randomIsing_2D(unsigned int len, unsigned int wid)
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
orderedIsing_2D(unsigned int len, unsigned int wid)
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
Ising_2D::Ising_2D(int size, double T, double J, double h) : 
len(size), wid(size), nSpins(len * wid), Ising_Base(T, J, h)
{
	this->spins = randomIsing_2D(len, wid);
	this->E = this->computeEnergy();
	this->M = this->computeMag();
}


Ising_2D::Ising_2D(int len, int wid, double T, double J, double h) : 
len(len), wid(wid), nSpins(len * wid), Ising_Base(T, J, h)
{
	this->spins = randomIsing_2D(len, wid);
	this->E = this->computeEnergy();
	this->M = this->computeMag();
}


int
Ising_2D::update(){
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
Ising_2D::size() const
{
	return this->nSpins;
}


inline const std::vector<int>&
Ising_2D::getSpins() const
{
	return this->spins;
}


double
Ising_2D::computeEnergy()
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
Ising_2D::computeMag(){
	double acc = 0;
	for(unsigned i = 0; i < this->nSpins; i++)
	{
		acc += this->spins[i];
	}
	return acc;
}


const int&
Ising_2D::operator()(unsigned int i, unsigned int j)
{
	if((i > this->len - 1) || (j > this->wid - 1))
	{
		throw "|  function: Ising_2D::operator ()  |  file: Ising_2D.cpp  |  error: Index out of bounds |";
	}
	return this->spins[i * this->len + j];
}


inline bool
Ising_2D::operator==(const Ising_2D &other)
{
	if((this->len != other.len) || (this->wid != other.wid))
	{
		return false;
	}
	return this->spins == other.spins;
}