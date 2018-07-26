/*
 * lattice.cpp

 *
 *  Created on: Apr 8, 2017
 *      Author: johnstanco
 */

#include "../include/lattice.hpp"


arma::mat
lattice::recipLat()
{
  arma::mat	_recip(dim, dim);
  
  //3D case, constructing b1, b2, b3, from a1, a2, a3.
  if (dim == 3){
    _recip.col(0) = 2 * pi * cross(vectors.col(1), vectors.col(2)) / cell_vol;
    _recip.col(1) = 2 * pi * cross(vectors.col(2), vectors.col(0)) / cell_vol;
    _recip.col(2) = 2 * pi * cross(vectors.col(0), vectors.col(1)) / cell_vol;
  }
  //2D case, constructing b1, b2 from a1, a2.
  else if (dim == 2){
    _recip(0, 0) = -vectors(1, 1);
    _recip(0, 1) = vectors(1, 0);
    _recip(1, 0) = vectors(0, 1);
    _recip(1, 1) = -vectors(0, 0);
    _recip.col(0) *= 2 * pi / cell_vol;
    _recip.col(1) *= 2 * pi / cell_vol;
  } else {
    throw "function recipLat : Improper lattice Dimension. Lattice must be of dimension 2 or 3.";
  }
  return _recip;
}


double
lattice::cellVolume(){
  if(dim == 1){
    return norm(vectors.col(0));
  }
  else if(dim == 2){
    return dot(vectors.col(0), vectors.col(1));
  }
  if(dim == 3){
    return dot(vectors.col(1), cross(vectors.col(2), vectors.col(3)));
  }
  return 0;
}



lattice::lattice(arma::mat &lat) : vectors(lat), dim(lat.n_cols), cell_vol(cellVolume()), recip(recipLat()) {}



lattice::lattice(const lattice &other) : vectors(other.vectors), dim(other.dim), cell_vol(other.cell_vol), recip(other.recip) {}