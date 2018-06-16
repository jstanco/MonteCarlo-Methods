//created by John Stanco on 5/25/18

#include "help.hpp"

#ifndef state_base_h
#define state_base_h


class state_base{
protected:
	virtual bool accept(double);
public:
	virtual int update() = 0;
	virtual int printData(){return 1;}
};


inline bool
state_base::accept(double logProb){
	if(logProb >= 0) return 1;
	return flip(exp(logProb));

}


#endif /* state_base_h */