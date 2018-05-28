//created by John Stanco on 5/11/18

//Basic Helper Functions for Metropolis Hastings Monte-Carlo;

#include "../include/help.hpp"

double min(double x1, double x2){
	if (x1 > x2){
		return x2;
	}
	return x1;
}


double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


int flip(double x){
	if(0 > x || x > 1){
		throw "|  function: flip  |  file: testState.cpp  |  error:  input value 'x' must be between 0 and 1  |";
	}
	if(fRand(0, 1) < x){
		return 1;
	}
	return 0;
}