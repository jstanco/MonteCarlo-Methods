//created by John Stanco 6.20.18

#include <iostream>
#include "../include/help.hpp"

template<class T>
int printArr(const T* arr, size_t len){
	for(size_t i = 0; i < len; i++){
		std::cout << arr[i] << std::endl;
	}
	return 1;
}


class discrete_pdf{
	friend class discrete_cdf;
private:
	double *pdf;
	size_t len;
	
	friend void swap(discrete_pdf& first, discrete_pdf& second){
		std::swap(first.len, second.len);
		std::swap(first.pdf, second.pdf);
	}

public:

	discrete_pdf(size_t len = 0) : len(len), pdf(len? new double[len]() : nullptr){}

	discrete_pdf(const double _pdf[], size_t len) : pdf(new double(*_pdf)), len(len) {
		memcpy(pdf, _pdf, sizeof(double) * len);
	}

	//copy constructor
	discrete_pdf(const discrete_pdf& other) : pdf(new double(*other.pdf)), len(other.len) {
		memcpy(pdf, other.pdf, sizeof(double) * len);
	}

	~discrete_pdf(){delete pdf;}

	int print(){return printArr(pdf, len);}

	discrete_pdf& operator=(discrete_pdf other){
		swap(*this, other);
		return *this;
	}

	double operator[](size_t i) const{
		return len? pdf[i] : 0;
	}

};


class discrete_cdf{
private:

	double *cdf;
	size_t len;

	friend void swap(discrete_cdf& first, discrete_cdf& second){
		std::swap(first.len, second.len);
		std::swap(first.cdf, second.cdf);
	}

	size_t bisect(double u, size_t, size_t);

public:

	discrete_cdf(size_t len = 0) : len(len), cdf(len? new double[len]() : nullptr){}

	discrete_cdf(const double _cdf[], size_t len) : cdf(new double(*_cdf)), len(len) {
		memcpy(cdf, _cdf, sizeof(double) * len);
	}

	//copy constructor
	discrete_cdf(const discrete_cdf& other) : cdf(new double(*other.cdf)), len(other.len) {
		memcpy(cdf, other.cdf, sizeof(double) * len);
	}

	discrete_cdf(const discrete_pdf& pdf) : cdf(new double[pdf.len]()), len(pdf.len) {
		cdf[0] = pdf[0];
		for(size_t i = 1; i < len; i++){
			cdf[i] = cdf[i - 1] + pdf[i];
		}
		normalize();
	}

	~discrete_cdf(){delete cdf;}

	size_t sample();

	int normalize(){
		if(!len || cdf[len - 1] == 1 ) return 0;
		for(size_t i = 0; i < len; i++){
			cdf[i] /= cdf[len - 1];
		}
		return 1;
	}

	int print(){return printArr(cdf, len);}

	discrete_cdf& operator=(discrete_cdf other){
		swap(*this, other);
		return *this;
	}

	double operator[](size_t i) const{
		return len? cdf[i] : 0;
	}
};


size_t
discrete_cdf::bisect(double u, size_t i_0, size_t i_f){
	size_t sub_len = i_f - i_0;
	size_t midpt = i_0 + sub_len / 2;
	if(sub_len < 2) return midpt + 1;
	if(u < cdf[midpt]) return bisect(u, i_0, midpt);
	return bisect(u, midpt, i_f);
}


//pass in the pdf table to the cdf...
//randomly sort table, and then walk through it ?
//1 -> find neighborhood
//compute pdf based on this neighborhood
//sample the pdf using cdf
//Challenge -> identifying neighborhood


size_t
discrete_cdf::sample(){
	double u = fRand(0, 1);
	return bisect(u, 0, len);
}

void space(){
	std::cout << std::endl;
}

int main(){

	size_t len = 10;
	double c[len];

	for(size_t i = 0; i < len; i++){
		c[i] = (i + 1) * (double)1 / len;
	}

	//printArr(c, len);
	discrete_pdf pdf(c, len);
	pdf.print();
	space();
	discrete_cdf cdf(pdf);
	cdf.print();
	space();

	std::cout << cdf.sample() << std::endl;

	

	//cdf.print();
	

	return 0;
}