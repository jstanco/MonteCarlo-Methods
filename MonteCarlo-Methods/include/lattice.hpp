/*
 * lattice.hpp

 *
 *  Created on: Apr 8, 2017
 *      Author: johnstanco
 */

#include "world"
#include <sstream>
#include </usr/local/include/armadillo>

#ifndef lattice_hpp
#define lattice_hpp

//Currently supports any output of wannier90;
class Lattice{
private:
  const arma::mat   vectors;
  const arma::mat   recip;
  const arma::vec   mass;
  const arma::vec   sigma;
  const arma::vec   eps;
  const double      cell_vol;
  const uint        dim;
  const std::string name_;
  double cellVolume();
  arma::mat recipLat();
public:
  //Lattice();
  Lattice(arma::mat &lat);
  Lattice(const Lattice &other);
};

#endif /* lattice_hpp */
