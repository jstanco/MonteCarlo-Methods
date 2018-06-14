//created by John Stanco on 5/25/18

#include "help.hpp"

#ifndef State_Base_Bones_h
#define State_Base_Bones_h


class State_Base_Bones{
protected:
	virtual bool accept(double);
public:
	virtual int update() = 0;
	virtual int printData(){return 1;}
};


inline bool
State_Base_Bones::accept(double logProb){
	if(logProb >= 0) return 1;
	return flip(exp(logProb));

}


#endif /* State_Base_Bones_h */