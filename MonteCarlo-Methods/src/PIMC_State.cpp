//created by John Stanco 6.2.18

#include "../include/pimc.hpp"

double hbar = 1;
double alpha = hbar * hbar;
std::ofstream ofs;


std::vector<bool> zeros(uint n){
	std::vector<bool> arr(n);
	for(uint i = 0; i < arr.size(); i++){
		arr[i] = false;
	}
	return arr;
}


int
pimc_state::startMove(const int pInd[], uint pSize){
	for(uint i = 0; i < pSize; i++){
		isMoved[pInd[i]] = true;
	}
	return 1;
}


int
pimc_state::endMove(const int pInd[], uint pSize){
	for(uint i = 0; i < pSize; i++){
		isMoved[pInd[i]] = false;
	}
	return 1;
}


pimc_state::pimc_state(double T, uint n_slice, 
					const arma::vec& masses, 
					const arma::vec& sigma, 
					const arma::vec& eps, 
					const slice& paths_mat) : 
					T(T), 
					dim(paths_mat.n_rows),
					n_part(paths_mat.n_cols),
					n_slice(n_slice),
					l_max(floor(log2(n_slice - 1))),
					clip_size(pow(2, floor(log2(n_slice - 1))) + 1),
					t0_max(n_slice - pow(2, floor(log2(n_slice - 1))) - 1),
					therm_wl(calcWl(T, n_slice, masses)),
					sigma(sigma), 
					eps(calcEps(T, n_slice, eps)),
					paths(calcPaths(n_slice, paths_mat)) {}

pimc_state::builder&
pimc_state::builder::setTemp(double T){
	_T = T;
	return *this;
}


pimc_state::builder&
pimc_state::builder::setSlices(uint n_slice){
	_n_slice = n_slice;
	return *this;
}


pimc_state::builder&
pimc_state::builder::setMasses(const arma::vec& masses){
	_masses = masses;
	return *this;
}


pimc_state::builder&
pimc_state::builder::setSigma(const arma::vec& sigma){
	_sigma = sigma;
	return *this;
}


pimc_state::builder&
pimc_state::builder::setEps(const arma::vec& eps){
	_eps = eps;
	return *this;
}


pimc_state::builder&
pimc_state::builder::setPaths(const arma::mat& paths){
	_paths = paths;
	return *this;
}


pimc_state
pimc_state::builder::build(){
	return pimc_state(_T, _n_slice, _masses, _sigma, _eps, _paths);
}


arma::vec
pimc_state::calcWl(double T, uint n_slice, const arma::vec& m){
	arma::vec _therm_wl(m.size());
	for(uint i = 0; i < m.size(); i++){
		_therm_wl(i) = pow(2 * m(i) * k * T * n_slice / alpha, -.5);
	}
	return _therm_wl;
}	


arma::vec
pimc_state::calcEps(double T, uint n_slice, const arma::vec& E){
	return arma::vec(E / (k * T * n_slice));
}


int
pimc_state::setBead(uint p_index, uint t, const bead& r){
	paths.tube(p_index, t) = r;
	return 1;
}


path_storage
pimc_state::calcPaths(uint n_slice, const slice& init){
	path_storage _path(init.n_cols, n_slice, init.n_rows);
	for(uint i = 0; i < init.n_cols; i++){
		for(uint j = 0; j < n_slice; j++){
			_path.tube(i, j) = init.col(i);
		}
	}
	isMoved = zeros(_path.n_rows);
	return _path;
}


int
pimc_state::moveBead(uint p_index, uint t, const arma::vec& dr){
	paths.tube(p_index, t) += dr;
	return 1;
}


int
pimc_state::moveParticle(uint p_index, const arma::vec& dr){
	for(uint t = 0; t < n_slice; t++){
		moveBead(p_index, t, dr);
	}
	return 1;
}


double
pimc_state::potential(uint p_index, const bead& r0, const bead& r1){
	double r = norm(r1 - r0);
	if(r < therm_wl(p_index)){ 
		double A = pow(sigma(p_index) / r, 6);
		return A * (A - 1);
	}
	return 0;
}


double
pimc_state::kinetic(uint p_index, const bead& r1, const bead& r2){;
	return pow(norm(r1 - r2), 2);
}


double
pimc_state::potentialChange(uint p_index, uint t, const bead& r0, const bead& r1){
	double dS = 0;
	for(uint i = 0; i < n_part; i++){
		if(i != p_index && !(isMoved[i] && i < p_index)){
			dS -= potential(p_index, r0, paths.tube(i, t));
			dS += potential(p_index, r1, paths.tube(i, t));
		}
	}
	return dS * eps(p_index);
}


double
pimc_state::kineticChange(uint p_index, uint t, const bead& r0, const bead& r1){
	double dS = 0;
	dS -= kinetic(p_index, r0, paths.tube(p_index, (t + 1) % n_slice));
	dS -= kinetic(p_index, r0, paths.tube(p_index, (t + n_slice - 1) % n_slice));
	dS += kinetic(p_index, r1, paths.tube(p_index, (t + 1) % n_slice));
	dS += kinetic(p_index, r1, paths.tube(p_index, (t + n_slice - 1) % n_slice));
	return dS / (2 * pow(therm_wl(p_index), 2));
}


//can precompute n, s, t, dt to avoid multiplications


double
pimc_state::midPoint(const int pInd[], uint pSize, uint t, uint s){
	double dS = 0;
	for(uint i = 0; i < pSize; i++){
		arma::vec dr = therm_wl(pInd[i]) * sqrt(s) * normDistVec(dim);
		bead newBead = .5 * (paths.tube(pInd[i], t + s) + paths.tube(pInd[i], t - s));
		dS += potentialChange(pInd[i], t, paths.tube(pInd[i], t), newBead + dr);
		setBead(pInd[i], t, newBead + dr);
	}
	return dS;
}


int
pimc_state::rev_bisect(const int pInd[], uint pSize, uint t0, uint _l, uint _l_max){
	uint t, dt, s, n;

	for(uint l = 1; l <= _l; l++){
		s = pow(2, _l_max - l);
		n = pow(2, l - 1);
		for(int m = 0; m < n; m++){
			dt = (2 * m + 1) * s;
			t = t0 + dt;
			for(uint i = 0; i < pSize; i++){
				setBead(pInd[i], t, old_paths.tube(i, m));
			}
		}
	}
	return 1;
}



int
pimc_state::bisect(const int pInd[], uint pSize, uint t0, uint _l_max){
	int n, s, t, dt;
	double dS = 0;
	for(uint l = 1; l <= _l_max; l++){
		s = pow(2, _l_max - l);
		n = pow(2, l - 1);
		for(int m = 0; m < n; m++){
			dt = (2 * m + 1) * s;
			t = t0 + dt;
			for(uint i = 0; i < pSize; i++){
				//std::cout << pInd[i] << "  " << m << std::endl;
				old_paths.tube(i, m) = paths.tube(pInd[i], t);
			}
			dS += midPoint(pInd, pSize, t, s);
		}
		if(!accept(-dS * s)){
			rev_bisect(pInd, pSize, t0, l, _l_max);
			return 0;
		}
	}
	return 1;
}


int
pimc_state::inverseCycle(const int cycle[], uint pSize, uint t){
	bead tmp = paths.tube(cycle[pSize - 1], t);
	for(uint i = pSize - 1; i > 0; i--){
		setBead(cycle[i], t, paths.tube(cycle[i - 1], t));
	}
	return setBead(cycle[0], t, tmp);
}



int
pimc_state::cycleEndpts(const int cycle[], uint pSize, uint t){
	bead tmp = paths.tube(cycle[0], t);
	for(uint i = 0; i < pSize - 1; i++){
		setBead(cycle[i], t, paths.tube(cycle[i + 1], t));
	}
	return setBead(cycle[pSize - 1], t, tmp);
}


int
pimc_state::permute(const int pInd[], uint pSize){

	int t0 = randInt(t0_max);
	int tf = t0 + clip_size - 1;
	old_paths = path_storage(pSize, clip_size, dim);
	cycleEndpts(pInd, pSize, tf);
	if(bisect(pInd, pSize, t0, l_max)){
		return 1;
	} 

	inverseCycle(pInd, pSize, tf);
	return 0;
}


//store permutation as cycle (k1 k2 k3 ... kn)
//This is pInd ... tells us what to do with each one
//We move from pInd[i] to pInd[i-1]

//should define the single-slice sampling to be of the form
//in which we sample from a gaussian first...



int
pimc_state::permutation(){
	//handle permutation logic
	int pSize = 4;
	int pInd[pSize];
	pInd[0] = 1;
	pInd[1] = 2;
	pInd[2] = 3;
	pInd[4] = 0;
	startMove(pInd, pSize);

	int rVal = permute(pInd, pSize);
	
	endMove(pInd, pSize);
	return rVal;
}

//var buildPermutation = function(particles)
//if(empty(particles)) return []
//i = uniformDraw(particles)
//if(good(i)) return [i].concat(buildPermutation(particles.remove[i]))
//return buildPermutation(particles)


int
pimc_state::bisection(){
	//want to ensure that the number of moves is 
	int moveSize = 1;//randInt(n_part);
	int pInd[moveSize];

	pInd[0] = randInt(n_part);

	//can create a random vector of particle moves
	int t0 = randInt(t0_max);

	startMove(pInd, moveSize);
	old_paths = path_storage(moveSize, clip_size, dim);
	
	int rVal = bisect(pInd, moveSize, t0, l_max);

	endMove(pInd, moveSize);
	return rVal;
}


double
pimc_state::checkParticleMove(uint p_index, const arma::vec& dr){
	double dS = 0;
	for(unsigned i = 0; i < n_slice; i++){
		bead tmp = paths.tube(p_index, i);
		dS += potentialChange(p_index, i, tmp, tmp + dr);
	}	
	return dS;
}


int
pimc_state::centerOfMass(){
	int p_index = randInt(n_part);
	arma::vec dr = uniDistVec(dim) * therm_wl(p_index);	
	double dS = checkParticleMove(p_index, dr);
	//printf("Change in action from moving particle %u:\t%f\n", p_index, dS);
	if(accept(-dS)){
		//printf("Particle move accepted!\n");
		moveParticle(p_index, dr);
		return 1;
	} else {
		//printf("Particle move rejected!\n");
	}
	return 0;
}


double
pimc_state::checkSingleBeadMove(uint p_index, uint t, const bead& dr){
	double dS = 0;
	arma::vec tmp = paths.tube(p_index, t);
	dS += potentialChange(p_index, t, tmp, tmp + dr);
	dS += kineticChange(p_index, t, tmp, tmp + dr);
	//printf("dS: %f\n", dS);
	return dS;
}


int
pimc_state::singleSlice(){
	int p_index = randInt(n_part);
	int t = randInt(1, n_slice - 1);
	arma::vec dr = uniDistVec(dim) * therm_wl(p_index);
	double dS = checkSingleBeadMove(p_index, t, dr);
	//printf("Change in action from moving bead %u, %u:\t%f\n", p_index, t, dS);
	if(accept(-dS)){
		//printf("Single-slice move accepted!\n");
		moveBead(p_index, t, dr);
		//update observables here
		return 1;
	} else {
		//printf("Single-slice move rejected!\n");
	}
	return 0;
}


//--------------------------------------------------------------//
//						   Public Methods						//
//--------------------------------------------------------------//


int
pimc_state::update(){
	//printf("--------------------------------------------\n");
	if(flip(.3)){
		return permutation();
	} else if (flip(.3)){
		return centerOfMass();
	} else {
		return bisection();	
	}
	//printf("--------------------------------------------\n");
	return 1;
}


arma::vec
pimc_state::getPath(uint p_index){
	//computes center of mass of particle (costly)
	arma::vec com(dim);
	com.zeros();
	for(uint i = 0; i < n_slice; i++){
		com += paths.tube(p_index, i);
	}
	return com / n_slice;
}


double
pimc_state::variance(uint p_index){
	arma::vec mean = getPath(p_index);
	double var = 0;
	for(uint i = 0; i < n_slice; i++){
		arma::vec curr = paths.tube(p_index, i);
		var += pow(norm(mean - curr), 2);
	}
	return var / n_slice;
}

int 
pimc_state::printToFile(std::string filename){
	ofs.open(filename, std::ios_base::app);
	for(uint i = 0; i < n_part; i++){
		for(uint j = 0; j < n_slice; j++){
			for(uint l = 0; l < dim; l++){
				ofs << paths(i, j, l) << "\t\t";
			}
			ofs << '\n';
		}
	}
	return 1;
}