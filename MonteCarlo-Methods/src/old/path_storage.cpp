//created by John Stanco 6.15.18

#include "../include/pimc.hpp"


path_storage::path_storage(uint n_slice, slice init) : paths(arma::cube(init.n_cols, n_slice, init.n_rows))
{
	for(uint i = 0; i < init.n_cols; i++){
		for(uint j = 0; j < n_slice; j++){
			paths.tube(i, j) = init.col(i);
		}
	}
}


path_storage::path_storage(uint n_part, uint n_slice, uint dim) : paths(arma::cube(n_part, n_slice, dim)) {}


int
path_storage::setBead(uint p_index, uint t, const bead& r){
	this->bead_at(p_index, t) = r;
	return 1;
}


int
path_storage::moveParticle(uint p_index, const arma::vec& dr){
	for(uint t = 0; t < paths.n_cols; t++){
		moveBead(p_index, t, dr);
	}
	return 1;
}


int
path_storage::moveBead(uint p_index, uint t, const arma::vec& dr){
	this->bead_at(p_index, t) += dr;
	return 1;
}


int
path_storage::replaceSlice(const int pInd[], uint pSize, const slice& newSlice){
	return 1;
}


bead&
path_storage::bead_at(uint p_index, uint t){
	return paths.tube(p_index, t);
}


path&
path_storage::path_at(uint p_index){
	return paths.row(p_index);
}


slice&
path_storage::slice_at(uint t){
	return paths.col(t);
}
