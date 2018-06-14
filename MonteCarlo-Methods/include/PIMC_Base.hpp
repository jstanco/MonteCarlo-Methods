//created by John Stanco 6.1.18


#ifndef PIMC_Base_h
#define PIMC_Base_h


#include "State_Base_Bones.hpp"


class PIMC_State : public State_Base_Bones{
protected:


public:
	PIMC_Base(){}
	PIMC_Base(double, uint, uint, uint);

	const double temp() 			const 	{return this->_T;}
	const uint n_slices()			const 	{return this->n_slice;}
	const uint n_particles()		const 	{return this->n_part;}
};


inline



#endif /* PIMC_Base_h */