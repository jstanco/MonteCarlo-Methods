//created by John Stanco on 5/10/18

#include "../include/testState.hpp"


//------------------------------------------------------//

double p(double x){
	return exp(-x);
}

//--------------------------------------------------------//


testState::testState(double x) : data(x) , prob(p(x)) {}


testState
testState::newState(){
	//generate new state based on some distribution.
	return *this;
}


double testState::priorProb(testState proposed){
	return 1;
}


double testState::targetProb(){
	return prob;
}


double
testState::acceptanceProb(testState proposed){
	double n = proposed.targetProb() * proposed.priorProb(*this);
	double d = this->prob * this->priorProb(proposed);
	if(d == 0){
		return 1;
	}
	return min(1, n / d);
}


double
testState::logProb(){
	return log(prob);
}


int
testState::accept(testState proposed){
	return (flip(this->acceptanceProb(proposed)));
}


testState*
testState::transition(){
	testState proposed = this->newState();
	if(accept(proposed)){
		return new testState(proposed);
	}
	return new testState(*this); //returns copy of original state? test to see if references work instead -> if I can also store references?
}