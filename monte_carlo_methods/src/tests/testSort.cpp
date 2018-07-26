//created by John Stanco on 5/18/18

#include <iostream>
#include <ctime>
#include <vector>

//testing sort


std::vector<int>
randomIntVector(unsigned int len){
	std::vector<int> ints(len);
	for(unsigned i = 0; i < len; i++){
		ints[i] = 2 * (rand() % 2) - 1;
	}
	return ints;
}


int
compare(int a, int b){
	return (a >= b);
}


int
greater(std::vector<int> arr1, std::vector<int> arr2){
	if(arr1.size() != arr2.size()){
		throw "|  function: greater  |  file: testSort.cpp  |  error:  input arrays must be of same length  |";
	}
	int len = arr1.size();
	for(unsigned i = 0; i < len; i++){
		if(arr1[i] != arr2[i]){
			return arr1[i] > arr2[i];
		}
	}
	return false;
}


//sort arrays by values inside array --> 


std::vector<int>
sortIntVector(std::vector<int> ints){
	std::sort(ints.begin(), ints.end(), compare);
	return ints;
}


template <class T>
int print(std::vector<T> arr1, std::vector<T> arr2, unsigned int len){
	for(unsigned i = 0; i < len; i++){
		std::cout << arr1[i] << "          " << arr2[i] << std::endl;
	}
	return 1;
}


int testArrayComparison(){
	srand(time(NULL));
	int len = 10;
	std::vector<int> spins1 = randomIntVector(len);
	std::vector<int> spins2 = randomIntVector(len);
	print(spins1, spins2, len);
	return greater(spins1, spins2);
}


int main(){
	printf("%d\n", testArrayComparison());
}