//created by John Stanco on 5/11/18

#include <random>
#include "mc_methods.hpp"

#ifndef HELP_H
#define HELP_H

int flip(double x = .5);
int rand_int(size_t, size_t);
int rand_int(size_t);
double min(double, double);
double fRand(double, double);
double plusMinus(double);
double randNorm();
arma::vec normDistVec(size_t);
arma::vec	uniDistVec(size_t);
arma::vec uniDistVec(const arma::vec&, const arma::vec&);


#endif /* HELP_H */