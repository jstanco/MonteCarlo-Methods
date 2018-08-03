# MonteCarlo-Methods

This is a software project designed to provide both functionality in solving problems involving magnetism, as well as a more general-purpose c++ API for other problems that might require a Markov-Chain monte-carlo approach.  (For example: a rational-rules model for learning in computational cognition, or a stochastic process in finance)

For the moment, The API consists of some templated MCMC classes for updating a user-defined state, as well as the namespace Signal, which provides a class for storing time-series data (real or complex) and manipulating it with the very fast FFTW3 Fourier Transform library.  The Signal namespace also includes some methods for statistical anlaysis of the data over time (mean, variance, cross-correlation, convolution, ...) with more soon to come.  

As a proof of functionality of the API, there is an implementation written for an n-demensinal square ising lattice.  This involves writing a class which functions as the 'State' which the MCMC object will then run.  The 'State' is a template class that would require certain methods, such as 'update(bool)' or 'prepare(int)', which allow the 'updater' to interface with the state it is updating.  

The Ising State consists of a general rank-n tensor object that stores the value of the spin at each lattice site (implemented as single dimensional vector length N1xN2x...xNd) implemented using the Armadillo API for linear algebra (http://arma.sourceforge.net/), along with a series of different child-classes that implement custom update procedures in accordance with different algorithms.  Currently, there are both single update (Metropolis-Hastings, Heat-Bath, Glauber) and clustering (Wolff) algorithms.  Note, while the heat-bath and Glauber procedures are different in theory, they do reduce to the same acceptance probability in the case of the Ising System.  I have included both for the sake of generality.  

I am in the process of the adding the worm algorithm.  I am also adapting the model to quantum spin systems, with the goal of simulating not just a square lattice, but triangular and polycrystalline (using a delaunay triangulation).  

The two basic procedures available are:
  1) MCMC run with the algorithm of choice: From the command line specify: lattice size, iterations, temperature, z-field and coupling constant J in [+1,-1].  This prints the expectation value of energy and abs(magnetization), along with time taken.
  
  2) From within main file (ising_temp_dependence.cpp), specify algorithm as well as range of temperatures to examine: energy, abs(mag), specific heat, and mag susceptibility.  These results currently agree nicely in 2D with Onsager's analytic solution, as well as with other MCMC simulations.
  
There is much more room to grow, and as mentioned before, all obervables from the run are stored as a time-series signal, and thus the analysis of the autocorellation times is a single function call away.  I will include an exponential fit function in the future.

Finally, one finds within the source code the skeleton for the next iteration of this project, which is a general purpose quantum PIMC (Path Integral Monte Carlo) simulation, designed to model molecular dynamics in either a fluid or equilibrium lattice configuration.  There was a working version of this about a month prior, although it is since deprecated as it did not utilize the data structure necessary to permute particles in a resonable amount of time.  Thus, it is being revamped as an extension of the Ising model.  

TO COMPILE:

The user must link their own BLAS,LAPACK,FFTW3,armadillo libraries.  The file "mc_methods.hpp" has paths to the approriate header files defined as macros.  One must substitute the paths to the headers from their own system.  The same must be done for the paths to the libraries defined in the makefile.  (This will all be made more streamlined in the future, but for the moment it is functional)

Once the paths for the headers are included, from within the monte_carlo_methods directory type in bash:
  
  cd src && make all
  
The binaries for the two routines are respectively: /bin/ising and /bin/ising_t_dependence.  
For more details on the ising model see: https://arxiv.org/pdf/1404.0209.pdf.  
As for the required libraries: 

  BLAS : Basic Linear Algebra Subroutines | https://www.openblas.net/ | There are MANY distributions/versions of this, but    OpenBlas is both open source and one of the fastest, using the GoToBlas source, which is used in supercomputing clusters.
  
  LAPACK : Linear Algebra PACKage | http://www.netlib.org/lapack/ | Couple of things: 1) This is not needed if OpenBlas is installed from the given website, as OpenBlas uses some LAPACK routines and thus includes its own version. 2) Note that the original source is written in FORTRAN77, and is usually linked to C++ via the 'extern "C"' keyword.  However, the main LAPACK implementations (including the one that comes with OpenBlas) include a separate implementation in C called LAPACKE.  Installing LAPACK will also install LAPACKE by default, although both the fortran and C libraries must be linked.
  
  Armadillo : | C++ library for linear algebra & scientific computing | (http://arma.sourceforge.net/) | This one also uses BLAS and LAPACK, and as such, will include its own versions of both.  These are not as optimized as OpenBlas, and thus one can go through the process of linking OpenBlas and its LAPACK(E) implementation to Armadillo if one sees fit.  
  
  FFTW3 : Fastest Fourier Transform (in the) West | This is separate from all of the others, and is used for ridiculously fast fourier transforms of both single and multidimenional data of both real and complex numbers.  It uses a type FFTW_complex which is really just an aliased double[2].  There is a way to link it such that the default FFTW_complex type is double \_Complex from the <complex.h> c library, although it is not assumed in this code.
  



