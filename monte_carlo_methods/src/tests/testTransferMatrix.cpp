//created by John Stanco on 10/12/18

//testing if transfer matrix agrees with analytic solution in thermodynamic limit.


#include "/usr/local/include/armadillo"


#define B 1
#define h .7
#define J 1


arma::mat
IsingTransferMatrix(){
	arma::mat T(2, 2);
	T 	<< exp(B * (J + h)) << exp(-B * (J + h)) << arma::endr
		<< exp(B * (h - J)) << exp(B * (J - h)) << arma::endr;
	return T;
}


arma::mat
myPow(arma::mat A, unsigned int n){
	arma::mat tmp = arma::eye(A.n_rows, A.n_cols);
	if(n > 0){
		for(unsigned i = 0; i < n; i++){
			tmp = tmp * A;
		}
		return tmp;
	}
	return tmp;
}


double 
IsingPartitionFunction(unsigned int L){
	arma::mat T_L = myPow(IsingTransferMatrix(), L);
	int z = 0;
	for(unsigned i = 0; i < 2; i++){
		for(unsigned j = 0; j < 2; j++){
			z += T_L(i, j);
		}
	}
	return z;
}


double IsingFreeEnergy(int L){
	return -(1 / B) * log(IsingPartitionFunction(L)) / L;
}




double analyticFreeEnergy(){
	return -(1 / B) * log( exp(B * J) * cosh(B * h) + sqrt(exp(2 * B * J) * pow(sinh(B * h), 2) + exp(-2 * B * J) ) );
}


double sinh(double x){
	return (exp(x) - exp(-x)) / 2;
}

double cosh(double x){
	return (exp(x) + exp(-x)) / 2;
}


arma::mat
testMatrixPow(){
	arma::mat T(2, 2);
	T 	<< exp(B * (J + h)) << exp(-B * (J + h)) << arma::endr
		<< exp(B * (h - J)) << exp(B * (J - h)) << arma::endr;
	arma::mat T1 = myPow(T, 5);
	arma::mat T2 = arma::pow(T, 2);
	return T1 - T * T * T * T * T;
}

int main(){
	std::cout << IsingFreeEnergy(1000) << std::endl;
	std::cout << analyticFreeEnergy() << std::endl;
	
	return 1;
}