//created by John Stanco on 5/10/18
//abstract base class used to implement state

#include <string>
#include "help.hpp"


#ifndef state_h
#define state_h

class state{
protected:
	double prob;
public:
	virtual double logProb() = 0;
	virtual state& transition() = 0;
	//virtual std::string printData() = 0;
	//virtual double priorProb(state&) = 0;
	//virtual double targetProb();
	//virtual double acceptanceProb(state& other;
	//virtual int accept(state& other);
};

#endif /* state_h */

// define quantum-mechanical system using hopping energies as prior probability distribution
// once we implement the structure of quantum mechanics, we can then locate larger-scale ground state of
// many-body systems?  Is this necessary?