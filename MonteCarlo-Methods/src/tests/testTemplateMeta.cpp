//created by John Stanco 6.14.18

#include <iostream>


//iterative
long long int fact_iter(uint n){
	long long int f = n;
	for(uint i = n - 1; i > 1; i--){
		f *= i;
	}
	return f;
}



//recursive
long long int fact_recur(uint n){
	if(n < 2) return 1;
	return n * fact_recur(n - 1);
}



//template
template<long long int n>
class template_fact{
public:
	enum { value = n * template_fact<n - 1>::value };
};

template<>
class template_fact<0>{ public: enum { value = 1 }; };



int printTime(std::string f, uint n, long long int fn, float t){
	std::cout 	<< f << " of " << n 
				<< " returned " << fn << ", took " 
				<< t / CLOCKS_PER_SEC << " seconds." << std::endl;
	return 1;
}


long long int runFact(uint n, uint iter, long long int(*f)(uint)){
	long long int fn = 1;
	for(uint i = 0; i < iter; i++){
		fn = f(n);
	}
	return fn;
}


int testFactTime(uint n, uint iter){
	long long int fn;
	clock_t t = clock();
	fn = runFact(n, iter, fact_iter);
	t = clock() - t;

	printTime("Iterative fact", n, fn, t);
	t = clock();
	fn = runFact(n, iter, fact_recur);
	t = clock() - t;

	printTime("Recursive fact", n, fn, t);
	return 1;
}



int main(){
	testFactTime(20, 1e7);
	clock_t t = clock();
	long long int fn;
	for(uint i = 0; i < 1e7; i++) {
		fn = template_fact<20>::value;
	}
	t = clock() - t;
	printTime("template fact", 20, fn, t);
	return 1;
}