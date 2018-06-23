//created by John Stanco on 5/11/18

#include <cmath>
#include <iostream>
#include <random>
#include "world"
#include "/usr/local/include/armadillo"

#ifndef help_h
#define help_h


int flip(double x = .5);
int randInt(uint a, uint b);
int randInt(uint a);
double min(double x1, double x2);
double fRand(double fMin, double fMax);
double plusMinus(double x);
double randNorm();
arma::vec normDistVec(int n);
arma::vec uniDistVec(int n);


#endif /* help_h */