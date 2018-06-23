//created by John Stanco 6.1.18


#ifndef pimc_h
#define pimc_h


#include "state_base.hpp"
#include "lattice.hpp"

typedef arma::cube path_storage;
typedef arma::mat path;
typedef arma::mat slice;
typedef arma::vec bead;


class pimc_state : public state_base{
protected:

	const uint 	 	l_max;
	const uint 		clip_size;
	const uint		t0_max;
	std::vector<bool> isMoved;
	path_storage 	paths;
	path_storage	old_paths;
	const arma::vec therm_wl;
	//parameters for Lennard-Jones potential (atom-dependent)
	const arma::vec sigma;
	const arma::vec eps;
	

	arma::vec calcWl(double, uint, const arma::vec&);
	arma::vec calcEps(double, uint, const arma::vec&);
	path_storage calcPaths(uint, const arma::mat&);

	virtual double potential(uint, const bead&, const bead&);
	virtual double kinetic(uint, const bead&, const bead&);

	virtual double potentialChange(uint, uint, const bead&, const bead&);
	virtual double kineticChange(uint, uint, const bead&, const bead&);

	virtual int setBead(uint, uint, const bead&);
	virtual int moveParticle(uint, const bead&);
	virtual int moveBead(uint, uint, const bead&);

	virtual double checkParticleMove(uint, const bead&);
	virtual double checkSingleBeadMove(uint, uint, const bead&);
	virtual double midPoint(const int[], uint, uint, uint);
	virtual int	cycleEndpts(const int[], uint, uint);
	virtual int inverseCycle(const int[], uint, uint);

	
	virtual int permute(const int[], uint);
	virtual int rev_bisect(const int[], uint, uint, uint, uint);
	virtual int	bisect(const int[], uint, uint, uint);

	int startMove(const int[], uint);
	int endMove(const int[], uint);

	virtual int permutation();
	virtual int bisection();
	virtual int	centerOfMass();
	virtual int singleSlice();

	pimc_state(double, uint, const arma::vec&, const arma::vec&, const arma::vec&, const arma::mat&);

public:

	const double T;
	const uint 	 n_slice;
	const uint 	 n_part;
	const uint 	 dim;

	class builder{
		private:
			double 		_T;
			uint 		_n_slice;
			arma::vec 	_masses;
			arma::vec 	_sigma;
			arma::vec 	_eps;
			arma::mat 	_paths;
		public:	
			builder(){}
			builder& setTemp(double);
			builder& setSlices(uint);
			builder& setMasses(const arma::vec&);
			builder& setSigma(const arma::vec&);
			builder& setEps(const arma::vec&);
			builder& setPaths(const arma::mat&);
			pimc_state build();//calls private constructor
	};

	virtual int update();
	virtual arma::vec getPath(uint);
	virtual double variance(uint);
	virtual int printToFile(std::string);

};


#endif /* pimc_h */