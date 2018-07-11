//created by John Stanco on 5/11/18

#include <cmath>
#include <iostream>
#include <random>
#include "world"
#include "/usr/local/include/armadillo"

#ifndef help_h
#define help_h

int 		flip( 		 double x = .5 	);
int 		rand_int( 	 size_t, size_t );
int 		rand_int( 	 size_t 		);
double 		min( 		 double, double );
double 		fRand( 		 double, double );
double 		plusMinus( 	 double 		);
double 		randNorm();
arma::vec 	normDistVec( size_t 		);
arma::vec 	uniDistVec(  size_t 		);


#endif /* help_h */