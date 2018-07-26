//created by John Stanco 7.24.18

#include "../../include/signal.hpp"
#include "../../include/timer.hpp"


cx_double cexp(double x){
	return cx_double(cos(x),sin(x));
}

cx_double cx_box(double x){
	if(x<-1 || x>1) { return cx_double(0,0); }
	return cx_double(1,0);
}


cx_double cx_lin(double x){
	if(x<-1 || x>1) { return cx_double(0,0); }
	return cx_double(1-x,0);
}


double box(double x){
	if(x<-1 || x>1) { return 0; }
	return 1;
}

double lin(double x){
	if(x<-1 || x>1) { return 0; }
	return 1-x;
}


void test_signal(){
	int N = pow(2, 8);
	double tmin,tmax;
	tmin=-6;
	tmax=6;

	Signal::signal f(&lin,tmin,tmax,N);
	Signal::signal g(&box,tmin,tmax,N);
	Signal::cx_signal h(&cx_lin,tmin,tmax,N);

	timer t;

	Signal::signal ff = autocorr(g);
	Signal::signal ff_2 = autocorr_slow(g);
	print(fftshift(ff),"../data/signal_conv.dat");
	print(fftshift(ff_2),"../data/signal_conv_slow.dat");
}


int main(){
	test_signal();
	return 1;
}