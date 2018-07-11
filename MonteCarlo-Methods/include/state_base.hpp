//created by John Stanco on 5/25/18

#include "help.hpp"

#ifndef state_base_h
#define state_base_h


class state_base{
protected:
public:
	virtual int update( bool burn = false ) = 0;
	virtual int print_data(){ return 1; }
};


#endif /* state_base_h */