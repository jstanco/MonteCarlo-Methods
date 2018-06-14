//created by John Stanco 6.2.18

#include "../include/PIMC_State.hpp"

double hbar = 1;
double alpha = hbar * hbar;
std::ofstream ofs;

//test if array works better...
static arma::vec NULL_VEC;


//PIMC_State::PIMC_State(double T, const arma::vec& masses, const arma::cube& pos) : 
//PIMC_Base(T, pos.n_rows, pos.n_cols, pos.n_slices), masses(masses), pos(pos) {}





PIMC_State::PIMC_State(double T, uint n_slice, 
					const arma::vec& masses, 
					const arma::vec& sigma, 
					const arma::vec& eps, 
					const arma::mat& pos_mat) : 
					T(T), 
					dim(pos_mat.n_rows),
					n_part(pos_mat.n_cols),
					n_slice(n_slice),
					l_max(floor(log2(n_slice - 1))),
					clip_size(pow(2, floor(log2(n_slice - 1))) + 1),
					t0_max(n_slice - pow(2, floor(log2(n_slice - 1))) - 1),
					masses(calcMasses(T, n_slice, masses)),
					sigma(sigma), 
					eps(calcEps(T, n_slice, eps)),
					pos(calcPos(n_slice, pos_mat)) {}

PIMC_State::builder&
PIMC_State::builder::setTemp(double T){
	_T = T;
	return *this;
}


PIMC_State::builder&
PIMC_State::builder::setSlices(uint n_slice){
	_n_slice = n_slice;
	return *this;
}


PIMC_State::builder&
PIMC_State::builder::setMasses(const arma::vec& masses){
	_masses = masses;
	return *this;
}


PIMC_State::builder&
PIMC_State::builder::setSigma(const arma::vec& sigma){
	_sigma = sigma;
	return *this;
}


PIMC_State::builder&
PIMC_State::builder::setEps(const arma::vec& eps){
	_eps = eps;
	return *this;
}


PIMC_State::builder&
PIMC_State::builder::setPos(const arma::mat& pos){
	_pos = pos;
	return *this;
}


PIMC_State
PIMC_State::builder::build(){
	return PIMC_State(_T, _n_slice, _masses, _sigma, _eps, _pos);
}


arma::vec
PIMC_State::calcMasses(double T, uint n_slice, const arma::vec& m){
	return arma::vec(m * k * T * n_slice / (2 * alpha));
}	


arma::vec
PIMC_State::calcEps(double T, uint n_slice, const arma::vec& E){
	return arma::vec(E / (k * T * n_slice));
}


arma::cube
PIMC_State::calcPos(uint n_slice, const arma::mat& pos_mat){
	arma::cube _pos(pos_mat.n_cols, n_slice, pos_mat.n_rows);
	for(uint i = 0; i < pos_mat.n_cols; i++){
		for(uint j = 0; j < n_slice; j++){
			_pos.tube(i, j) = pos_mat.col(i);
		}
	}
	return _pos;
}


double
PIMC_State::potential(uint pIndex, const arma::vec& r0, const arma::vec& r1){
	double A = pow(sigma(pIndex) / norm(r1 - r0), 6);
	//std::cout << norm(r1 - r0) << std::endl;
	//std::cout << A << "\n" << std::endl;
	return A * (A - 1);
}


double
PIMC_State::kinetic(uint pIndex, const arma::vec& r1, const arma::vec& r2){;
	return pow(norm(r1 - r2), 2);
}


double
PIMC_State::potentialChange(uint pIndex, uint t, const arma::vec& r0, const arma::vec& r1){
	double dS = 0;
	for(uint i = 0; i < n_part; i++){
		if(i != pIndex){
			dS -= potential(pIndex, r0, pos.tube(i, t));
			dS += potential(pIndex, r1, pos.tube(i, t));
		}
	}
	return dS * eps(pIndex);
}


double
PIMC_State::kineticChange(uint pIndex, uint t, const arma::vec& r0, const arma::vec& r1){
	double dS = 0;
	dS -= kinetic(pIndex, r0, pos.tube(pIndex, (t + 1) % n_slice));
	dS -= kinetic(pIndex, r0, pos.tube(pIndex, (t + n_slice - 1) % n_slice));
	dS += kinetic(pIndex, r1, pos.tube(pIndex, (t + 1) % n_slice));
	dS += kinetic(pIndex, r1, pos.tube(pIndex, (t + n_slice - 1) % n_slice));
	return dS * masses(pIndex);
}


int
PIMC_State::replaceBead(uint pIndex, uint t, const arma::vec& r){
	pos.tube(pIndex, t) = r;
	return 1;
}


arma::vec
PIMC_State::midPoint(uint pIndex, uint t, uint stride){
	//std::cout << stride << std::endl;
	arma::vec r = .5 * (pos.tube(pIndex, t + stride) + pos.tube(pIndex, t - stride));
	return r + normDistVec(dim) * sqrt(((double)stride / 4) / masses(pIndex));
}


//particle swap
//Choose subpath of length t -> t + m
//Choose l in [2, N]
//compute permutation table


//walk through table
//perform permutation of endpoints, perform bisection procedure to test if newly generated 




/*
int
PIMC_State::permutation(){
	int _n_Part = randInt(n_part);
	int _l_max = randInt(l_max);
	int t0 = randInt(t0_max(_l_max));
	int tf = 
	//define t0_max on the fly?
	//can create table of all t_maxes for the different possible l_maxes
	//in general, log2()
}
*/




//can precompute n, s, t, dt to avoid multiplications
int
PIMC_State::bisection(){
	int pIndex = randInt(n_part);
	int n, s, t, dt; int t0 = randInt(t0_max);
	double dS = 0;
	bool accepted = true;
	arma::mat midPts(dim, clip_size);
	for(uint l = l_max; l > 0; l--){
		n = pow(2, l_max - l);
		s = pow(2, l - 1);
		for(int m = 0; m < n; m++){
			dt = (2 * m + 1) * s;
			t = t0 + dt;
			midPts.col(dt) = midPoint(pIndex, t, s);
			dS += potentialChange(pIndex, t, pos.tube(pIndex, t), midPts.col(dt));
		}
		//printf("Change in potential action from level %u: %f\n", l + 1, dS);
		if(!accept(-dS * s)){
			//printf("Bisection rejected on level: %u\n", l + 1);
			accepted = false;
			break;
		}
	}
	if(accepted){
		for(uint t = 1; t < clip_size - 1; t++){
			replaceBead(pIndex, t0 + t, midPts.col(t));
		}
		//update parameters
	}
	return 1;
}


double
PIMC_State::checkParticleMove(uint pIndex, const arma::vec& dr){
	double dS = 0;
	for(unsigned i = 0; i < n_slice; i++){
		arma::vec tmp = pos.tube(pIndex, i);
		dS += potentialChange(pIndex, i, tmp, tmp + dr);
	}	
	return dS;
}


int
PIMC_State::moveParticle(uint pIndex, const arma::vec& dr){
	for(unsigned i = 0; i < n_slice; i++){
		moveBead(pIndex, i, dr);
	}
	return 1;
}


int
PIMC_State::centerOfMass(){
	int pIndex = randInt(n_part);
	arma::vec dr = .5 * uniDistVec(dim) / sqrt(masses(pIndex));	
	double dS = checkParticleMove(pIndex, dr);
	//printf("Change in action from moving particle %u:\t%f\n", pIndex, dS);
	if(accept(-dS)){
		//printf("Particle move accepted!\n");
		moveParticle(pIndex, dr);
		//update observables
	} else {
		//printf("Particle move rejected!\n");
	}
	return 1;
}


double
PIMC_State::checkSingleBeadMove(uint pIndex, uint t, const arma::vec& dr){
	double dS = 0;
	arma::vec tmp = pos.tube(pIndex, t);
	dS += potentialChange(pIndex, t, tmp, tmp + dr);
	dS += kineticChange(pIndex, t, tmp, tmp + dr);
	//printf("dE: %f\n", dE);
	return dS;
}


int
PIMC_State::moveBead(uint pIndex, uint t, const arma::vec& dr){
	pos.tube(pIndex, t) += dr;
	return 1;
}


int
PIMC_State::singleSlice(){
	int pIndex = randInt(n_part);
	int t = randInt(1, n_slice - 1);
	arma::vec dr  = .5 * uniDistVec(dim) / sqrt(masses(pIndex));
	double dS = checkSingleBeadMove(pIndex, t, dr);
	//printf("Change in action from moving bead %u, %u:\t%f\n", pIndex, t, dS);
	if(accept(-dS)){
		//printf("Single-slice move accepted!\n");
		moveBead(pIndex, t, dr);
		//update observables here
	} else {
		//printf("Single-slice move rejected!\n");
	}
	return 1;
}

//TODO

//define a multislice /non-bisection move?
//get rid of the single-slice?
//use just bisection and COM?
//A cyclic permuation of many 

//--------------------------------------------------------------//
//						   Public Methods						//
//--------------------------------------------------------------//


int
PIMC_State::update(){
	//printf("--------------------------------------------\n");
	if(flip(.3)){
		bisection();
	} else if (flip(.3)){
		centerOfMass();
	} else {
		singleSlice();	
	}
	//printf("--------------------------------------------\n");
	return 1;
}


arma::vec
PIMC_State::getPos(uint pIndex){
	//computes center of mass of particle (costly)
	arma::vec com(dim);
	com.zeros();
	for(uint i = 0; i < n_slice; i++){
		com += pos.tube(pIndex, i);
	}
	return com / n_slice;
}


double
PIMC_State::variance(uint pIndex){
	arma::vec mean = getPos(pIndex);
	double var = 0;
	for(uint i = 0; i < n_slice; i++){
		arma::vec curr = pos.tube(pIndex, i);
		var += pow(norm(mean - curr), 2);
	}
	return var / n_slice;
}

int 
PIMC_State::printToFile(std::string filename){
	ofs.open(filename, std::ios_base::app);
	for(uint i = 0; i < n_part; i++){
		for(uint j = 0; j < n_slice; j++){
			for(uint l = 0; l < dim; l++){
				ofs << pos(i, j, l) << "\t\t";
			}
			ofs << '\n';
		}
	}
	return 1;
}