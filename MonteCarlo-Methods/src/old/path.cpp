//created by John Stanco 6.21.18

//potential data structure for worm algorithm

#include "../include/help.hpp"

class path_storage{

private:


public:
	arma::field<arma::vec> beads;

};


class bead{

private:

public:
	const size_t slice;
	const size_t p_index;
	bead* next;
	bead* prev;

	bead(size_t i = 0, size_t j = 0) : p_index(i), slice(j) {}
	int print(){std::cout << slice << "  " << p_index << std::endl; return 1;}
};

//have array of path...

class path{

public:
	size_t n_slice;
	bead* m; 	
	bool isWorm;


};







int main(){


	path p;
	bead *beads = new bead[5];

	p.m = beads;

	for(size_t i = 0; i < 5; i++){
		beads[i].next = &beads[(i + 1) % 5];
		beads[i].prev = &beads[(i - 1 + 5) % 5];
	}

	for(size_t i = 0; i < 5; i++){
		beads[i].next->print();
	}


	//contains array of paths -> no std::vector

	delete [] beads;

	return 1;
}


//class slice
//holds array of beads
//each bead should realistically know which bin it is in...

//in lattice formalism -> have structure lattice
//bin will contain pointer to the lattice
//


//re-implement using the 


//This is a data structure to properly implement the worm algorithm
//It is constructed in such a way

//In order to create these, we should have to realistically

//pimc_state<boltzmannon>
//pimc_state<fermion>






//Is it even necessary to create a class for path?

//could we simply have path_storage that stores all positions in a tensor,
//Tensor stores all of the next/prev info
//All that is really added is the addition of which 'bucket' each is stored in
//Bucket is simply represented by array...
//particles are initialized some way in the bucket
//How do we know that we are moving to a new bucket?
//Each bucket is characterized by a particular index?
//The corners define a bucket?
//When a move is made, how are we moving from one bucket to next?
//Define bucket modulo some characteristic length...

//Very similar to lattice, but particles are not repeated periodically, but rather
//can move freely between buckets
//When a particle is displaced, check how many times the particle as 

//should have two different modes : free and lattice
//paths provide extra layer of information?
//not needed, all adjacency info is stored in the particle table?
//This would be faster, as it would not ever copy vectors, but rather just copy a pointer!
//For permutations, we would not have to move the data, just swap pointers.
/*
	Is it possible to store all positions by pointers.
	Bucket is a pointer to the beads, handles which neighborhood they are in
	Instead, it would make more sense to just store a bucket index, and represent it virtually.
	Buckets store beads  -> beads store indices---> Bead gets moved to a different particle?
	Do the beads store the actual particle, or do they store the location for the particle?

	When a particle is swapped, its particle index is changed, along with its next variables
	When we update the system, sometimes the pointers get moved around, but no position vectors ever
	There are certain properties of each particle
	There are certain types of each particle...If particles are identical, can two different types of 
	particle be swapped?

	In that case, there would definitely be a change in the potential energy.
	Thus, it must only be particles of that same type that get get swapped

	should be a function 'is_swappable(bead&, bead&)'  this will tell us if the beads are
	of same type...  


	Store the beads in buckets.  The buckets contain a std::vector of beads.  
	From initialization, each bead contains a shallow copy prev/next element

	Buckets represent 

	//better idea:

	The entire structure uses the field class from arma.  At each location is a vector of particle indices
	The index of the bucket is its lattice coordinates
	So if we are using the free system, then the lattice is cubic and defined by some equidistant size
	However, if we are using a particular lattice, then there will be a lattice object that is used
	To define the coordinates and 'cell' of the lattice...Some modulo operation?
	Essentially, one can take the coordinates and determine the actual lattice coords, which would correspond
	to the bucket indices.

	*/
//


