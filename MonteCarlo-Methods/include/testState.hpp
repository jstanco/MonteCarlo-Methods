//created by John Stanco on 5/10/18

#include <cmath>
#include <iostream>
#include "state.hpp"

#ifndef testState_h
#define testState_h


class testState : public state{
protected:
	double data;
	double prob;
	//void calcTargetProb()
public:
	testState(double);
	testState newState();
	double priorProb(testState);
	double targetProb();
	double logProb();
	double acceptanceProb(testState);
	int accept(testState);
	testState* transition();
};


#endif /* testState_h */