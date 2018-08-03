//created by John Stanco on 5/25/18

#include "help.hpp"

#ifndef state_base_h
#define state_base_h

class state_base{
protected:
public:
	virtual void update(const bool burn = false) = 0;
	virtual void print_data() {}
};

#endif /* state_base_h */