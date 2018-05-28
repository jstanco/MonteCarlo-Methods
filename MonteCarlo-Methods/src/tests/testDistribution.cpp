//created by John Stanco on 5/19/18

//testing Distribution class

#include "../include/Distribution.hpp"
#include <iostream>
#include <ctime>
#include <cmath>


std::vector<int*>
randomIntVector(unsigned int len){
	std::vector<int*> ints(len);
	for(unsigned i = 0; i < len; i++){
		ints[i] = new int(rand() % 10);
		//ints[i] = 2 * (rand() % 2) - 1;
	}
	return ints;
}


template <class T>
int print(std::vector<T> arr1, std::vector<T> arr2, unsigned int len){
	for(unsigned i = 0; i < len; i++){
		std::cout << arr1[i] << "          " << arr2[i] << std::endl;
	}
	return 1;
}


template <class T>
int print(std::vector<T> arr1, unsigned int len){
	for(unsigned i = 0; i < len; i++){
		std::cout << arr1[i] << std::endl;
	}
	return 1;
}


template <class T>
int print(Distribution<T> d){
	unsigned int len = d.size();
	for(unsigned i = 0; i < len; i++){
		std::cout << d.data(i) << "       " << d.prob(i) << std::endl;
	}
	return 1;
}


template<class T>
double mean(T x){
	return x;
}


template<class T>
double secondMoment(T x){
	return x * x;
}


int testDistribution(){
	srand(time(NULL));
	int len = 10;
	std::vector<int*> tmp = randomIntVector(len);
	std::vector<int*> &ints = tmp;
	Distribution<int> d(ints);
	printf("Std. dev of distribution:  %f\n", sqrt(d.expVal(secondMoment) - pow(d.expVal(mean), 2) ) );
	print(d);
	return 1;
}


int main(){
	return testDistribution();
}