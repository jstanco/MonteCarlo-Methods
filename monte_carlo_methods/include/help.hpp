//created by John Stanco on 5/11/18

#include <random>
#include "mc_methods.hpp"

#ifndef HELP_H
#define HELP_H

int flip(const double x = .5);
int rand_int(const size_t, const size_t);
int rand_int(const size_t);
double min(const double, const double);
double fRand(const double, const double);
double plusMinus(const double);
double randNorm();
arma::vec normDistVec(const size_t);
arma::vec	uniDistVec(const size_t);
arma::vec uniDistVec(const arma::vec&, const arma::vec&);


#endif /* HELP_H */