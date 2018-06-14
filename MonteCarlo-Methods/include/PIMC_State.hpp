//created by John Stanco 6.1.18


#ifndef PIMC_State_h
#define PIMC_State_h


#include "State_Base_Bones.hpp"
#include "Lattice.hpp"


class PIMC_State : public State_Base_Bones{
protected:

	const uint 	 	l_max;
	const uint 		clip_size;
	const uint		t0_max;
	arma::cube 		pos;
	const arma::vec masses;
	//parameters for Lennard-Jones potential (atom-dependent)
	const arma::vec sigma;
	const arma::vec eps;
	

	arma::vec calcMasses(double, uint, const arma::vec&);
	arma::vec calcEps(double, uint, const arma::vec&);
	arma::cube calcPos(uint, const arma::mat&);

	virtual double potential(uint, const arma::vec&, const arma::vec&);
	virtual double kinetic(uint, const arma::vec&, const arma::vec&);

	virtual double potentialChange(uint, uint, const arma::vec&, const arma::vec&);
	virtual double kineticChange(uint, uint, const arma::vec&, const arma::vec&);

	virtual double checkParticleMove(uint, const arma::vec&);
	virtual double checkSingleBeadMove(uint, uint, const arma::vec&);
	virtual arma::vec midPoint(uint, uint, uint);

	virtual int replaceBead(uint, uint, const arma::vec&);
	virtual int moveParticle(uint, const arma::vec&);
	virtual int moveBead(uint, uint, const arma::vec&);

	virtual int bisection();
	virtual int	centerOfMass();
	virtual int singleSlice();

	PIMC_State(double, uint, const arma::vec&, const arma::vec&, const arma::vec&, const arma::mat&);

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
			arma::mat 	_pos;
		public:	
			builder(){}
			builder& setTemp(double);
			builder& setSlices(uint);
			builder& setMasses(const arma::vec&);
			builder& setSigma(const arma::vec&);
			builder& setEps(const arma::vec&);
			builder& setPos(const arma::mat&);
			PIMC_State build();//calls private constructor
	};

	PIMC_State();
	 //undefined rn

	virtual int update();
	virtual arma::vec getPos(uint);
	virtual double variance(uint);
	virtual int printToFile(std::string);

};


#endif /* PIMC_State_h */


//TODO -> make PIMC state one
//			Derive Fermionic and Bosonic systems from this
//			(update moves are roughly the same, just add permutations