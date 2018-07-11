//created by john stanco 6.18.18

#include "../include/help.hpp"


//Note, testing purposes only, no bounds-checking



int printArr(const int arr[], size_t len){
	for(size_t i = 0; i < len; i++){
		printf("%i\n", arr[i]);
	}
	return 1;
}



class cyclic_permutation{
	size_t len;
	int *cycle; 
public:
	cyclic_permutation(int *p, size_t len) : cycle(new int(*p)), len(len) {
		memcpy(cycle, p, sizeof(int) * len);
	}
	int print(){return printArr(cycle, len);}	
};



int
cycleEndpts(const int cycle[], size_t pSize, int num[]){
	int tmp = num[cycle[pSize - 1]];
	for(size_t i = pSize - 1; i > 0; i--){
		num[cycle[i]] = num[cycle[i - 1]];
	}
	num[cycle[0]] = tmp;
	return 1;
	
}



int
inverseCycle(const int cycle[], size_t pSize, int num[]){
	//print(num, 10);
	int tmp = num[cycle[0]];
	for(size_t i = 0; i < pSize - 1; i++){
		num[cycle[i]] = num[cycle[i + 1]];
	}
	num[cycle[pSize - 1]] = tmp;
	return 1;
}


int fillInOrder(int arr[], size_t len){
	for(size_t i = 0; i < len; i++){
		arr[i] = i;
	}
	return 1;
}


int swap(int arr[], uint i, uint j){
	int tmp = arr[i];
	arr[i] = arr[j];
	arr[j] = tmp;
	return 1;
}


int uniformDraw(const int arr[], size_t len){
	return arr[rand_int(len)];
}


int randomCycle(int cycle[], size_t pSize, size_t len){
	int pool[len];
	fillInOrder(pool, len);
	for(size_t i = 0; i < pSize; i++){
		int P_i = uniformDraw(pool, len - i);
		swap(pool, P_i, len - i - 1);
		cycle[i] = P_i;
	}
	return 1;
}


int testCycle1(){
	size_t len = 10;
	int num[len];
	fillInOrder(num, len);

	size_t pSize = 3;
	int cycle[pSize];
	randomCycle(cycle, pSize, len);
	//printf("%lx-cycle:\n", pSize);
	//printArr(cycle, pSize);
	return cycleEndpts(cycle, pSize, num);
}




int main(){
	size_t len = 4;
	int p[] = {1, 2, 3, 4};

	cyclic_permutation sigma(p, len);
	sigma.print();



	return 1;
}