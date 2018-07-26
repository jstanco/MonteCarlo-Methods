#ifndef MC_METHODS

/// All paths are specific to user
/// To link BLAS, LAPACK, FFTW, user must edit path both here and in /src/makefile

#ifndef LAPACKE_PATH
#define LAPACKE_PATH "/Users/johnstanco/PhysicsSoftware/lapack-3.8.0/LAPACKE/include/"
#define LAPACKE_HEADER "/Users/johnstanco/PhysicsSoftware/lapack-3.8.0/LAPACKE/include/lapacke.h"
#define LAPACKE_CONFIG_HEADER "/Users/johnstanco/PhysicsSoftware/lapack-3.8.0/LAPACKE/include/lapacke_config.h"
#define LAPACKE_MANGLING_HEADER "/Users/johnstanco/PhysicsSoftware/lapack-3.8.0/LAPACKE/include/lapacke_mangling.h"
#define LAPACKE_UTILS_HEADER "/Users/johnstanco/PhysicsSoftware/lapack-3.8.0/LAPACKE/include/lapacke_utils.h"
#endif /* LAPACKE_PATH */

#ifndef BLAS_PATH
#define BLAS_PATH "/Users/johnstanco/PhysicsSoftware/OpenBlas/"
#define BLAS_HEADER "/Users/johnstanco/PhysicsSoftware/OpenBlas/cblas.h"
#endif /* BLAS_PATH */

#ifndef FFTW_PATH
#define FFTW_PATH "/usr/local/include/"
#define FFTW_HEADER "/usr/local/include/fftw3.h"
#endif /* FFTW_PATH */

#ifndef ARMA_PATH
#define ARMA_PATH "/usr/local/include"
#define ARMA_HEADER "/usr/local/include/armadillo"
#endif /* ARMA_PATH */

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <complex.h>
#include <complex>

//#include LAPACKE_HEADER
//#include LAPACKE_CONFIG_HEADER
//#include LAPACKE_MANGLING_HEADER
//#include LAPACKE_UTILS_HEADER
#include BLAS_HEADER
#include FFTW_HEADER
#include ARMA_HEADER



extern double k;
extern double hBar;

#endif /* MC_METHODS */