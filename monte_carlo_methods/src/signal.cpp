//created by John Stanco 7.27.18

#include "../include/Signal.hpp"

/////////////////////// REAL Signal ////////////////////////

double Signal::mean(const TSignal<double>& f){
		double mu = 0;
		for(size_t i=0;i<f.len;i++){
			mu += f(i);
		}
		return mu / f.len;
	}


double Signal::var(const TSignal<double>& f){
	double sigma = 0;
	double mu = mean(f);
	for(size_t i=0;i<f.len;i++){
		sigma += pow(f(i) - mu,2);
	}
	return sigma / f.len;
}


void Signal::print(const Signal::TSignal<double>& f){
	for(size_t i=0;i<f.n_sample();i++){
		printf("t: %lf\tf: %lf\n", f.time(i),f(i));
	}
}


void Signal::print(const Signal::TSignal<double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for(size_t i=0;i<f.n_sample();i++){
			fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i)); }
	}
}

void Signal::print(const Signal::TSignal<double> &f, const std::string& filename){
	FILE *pFile = fopen(filename.c_str(), "w");
	if(pFile){ for(size_t i=0;i<f.n_sample();i++){
			fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i)); }
	}
}


/// DCT-II - input data must be real and symmetric
Signal::TSignal<double> Signal::fft(const Signal::TSignal<double>& f){
	int arrsize = f.len*sizeof(double);
	double *in = (double*)fftw_malloc(arrsize);
	double *out = (double*)fftw_malloc(arrsize);
	fftw_plan p = fftw_plan_r2r_1d(f.len, in, out, FFTW_REDFT10,FFTW_ESTIMATE);
	memcpy(in,f.data,arrsize);
	fftw_execute(p);
	TSignal<double> s(out,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return s;
}


Signal::TSignal<double> Signal::ifft(const Signal::TSignal<double>& f){
	int arrsize = f.len*sizeof(double);
	double *in = (double*)fftw_malloc(arrsize);
	double *out = (double*)fftw_malloc(arrsize);	//exploits fact that fft of real function is hermitian
	fftw_plan p = fftw_plan_r2r_1d(f.len, in, out,FFTW_REDFT01,FFTW_ESTIMATE);
	memcpy(in,f.data,arrsize);
	fftw_execute(p);
	TSignal<double> s(out,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return s;
}


/// IDCT-II - input data must be real and symmetric
Signal::TSignal<cx_double> Signal::fft_r2c(const Signal::TSignal<double>& f){
	int in_arrsize = f.len*sizeof(double);
	int out_arrsize = (f.len/2+1)*sizeof(fftw_complex);
	double *in = (double*)fftw_malloc(in_arrsize);
	fftw_complex *out = (fftw_complex*)fftw_malloc(out_arrsize);
	fftw_complex *rslt = (fftw_complex*)fftw_malloc(f.len*sizeof(fftw_complex));

	fftw_plan p = fftw_plan_dft_r2c_1d(f.len, in, out, FFTW_ESTIMATE);
	memcpy(in,f.data,in_arrsize);
	fftw_execute(p);

	/// FFTW only computes 1st half, use hermiticity for 2nd half | TODO: CHECK
	memcpy(rslt,out,out_arrsize);
	for(size_t i = f.len/2+1;i<f.len;i++){
		assign_conj(rslt[i],rslt[f.len-i]);
	}

	TSignal<cx_double> s(rslt,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(rslt);
	return s;
}


Signal::TSignal<double> Signal::ifft_c2r(const Signal::TSignal<cx_double>& f){
	int in_arrsize = (f.len/2+1)*sizeof(fftw_complex);
	int out_arrsize = f.len*sizeof(double);
	fftw_complex *in = (fftw_complex*)fftw_malloc(in_arrsize);
	double *out = (double*)fftw_malloc(out_arrsize);
	fftw_plan p = fftw_plan_dft_c2r_1d(f.len, in, out, FFTW_ESTIMATE);
	memcpy(in,f.data,in_arrsize);
	fftw_execute(p);
	TSignal<double> s(out,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return s;
}


Signal::TSignal<double> Signal::convolve(const Signal::TSignal<double>& l, const Signal::TSignal<double>& r){
	int n = l.len;
	int N = (n/2+1);
	int in_arrsize = n*sizeof(double);
	int out_arrsize = N*sizeof(fftw_complex);

	double *in = (double*)fftw_malloc(in_arrsize);
	fftw_complex *lout = (fftw_complex*)fftw_malloc(out_arrsize);
	fftw_complex *rout = (fftw_complex*)fftw_malloc(out_arrsize);
	fftw_complex *rslt = (fftw_complex*)fftw_malloc(out_arrsize);

	fftw_plan pl = fftw_plan_dft_r2c_1d(n, in, lout, FFTW_ESTIMATE);
	fftw_plan pr = fftw_plan_dft_r2c_1d(n, in, rout, FFTW_ESTIMATE);
	fftw_plan pc = fftw_plan_dft_c2r_1d(n, rslt, in, FFTW_ESTIMATE);

	memcpy(in,l.data,in_arrsize);
	fftw_execute(pl);
	memcpy(in,r.data,in_arrsize);
	fftw_execute(pr);

	mult_arr(lout,rout,rslt,N);

	fftw_execute(pc);
	TSignal<double> s(in,l.t,n,Signal::SAME);
	fftw_free(in);
	fftw_free(lout);
	fftw_free(rout);
	fftw_free(rslt);
	fftw_destroy_plan(pl);
	fftw_destroy_plan(pr);
	fftw_destroy_plan(pc);
	return s;
}


Signal::TSignal<double> Signal::crosscorr(const Signal::TSignal<double>& l, const Signal::TSignal<double>& r){
	return Signal::convolve(l,r);
}

Signal::TSignal<double> Signal::autocorr(const Signal::TSignal<double>& f){
	return Signal::crosscorr(f,f);
}


///////////////////// COMPLEX Signal ///////////////////////


void Signal::print(const Signal::TSignal<cx_double>& f){
	for(size_t i=0;i<f.n_sample();i++){
		printf("t: %lf\tre: %lf\tim: %lf\n", f.time(i),f(i).real(),f(i).imag());
	}
}

void Signal::print_real(const Signal::TSignal<cx_double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for(size_t i=0;i<f.n_sample();i++){
			fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i).real()); }
	}
}

void Signal::print_imag(const Signal::TSignal<cx_double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for(size_t i=0;i<f.n_sample();i++){
			fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i).imag()); }
	}
}


Signal::TSignal<cx_double> Signal::fft(const TSignal<cx_double>& f){
	int arrsize = sizeof(fftw_complex)*f.len;
	fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *out = (fftw_complex*)fftw_malloc(arrsize);
	fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	memcpy(in,f.data,arrsize);
	fftw_execute(p);
	TSignal<cx_double> s(out,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return s;
}


Signal::TSignal<cx_double> Signal::ifft(const TSignal<cx_double>& f){
	int arrsize = sizeof(fftw_complex)*f.len;
	fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *out = (fftw_complex*)fftw_malloc(arrsize);
	fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	memcpy(in,f.data,arrsize);
	fftw_execute(p);
	TSignal<cx_double> s(out,f.t,f.len,Signal::DUAL);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return s;
}


/// assume Signal::SAME length for now
Signal::TSignal<cx_double> Signal::convolve(const TSignal<cx_double>& l, const TSignal<cx_double>& r){
	int N = l.len;
	int arrsize = N*sizeof(fftw_complex);

	fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *lout = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *rout = (fftw_complex*)fftw_malloc(arrsize);

	fftw_plan pl = fftw_plan_dft_1d(N,in,lout,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan pr = fftw_plan_dft_1d(N,in,rout,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan pc = fftw_plan_dft_1d(N,in,lout,FFTW_BACKWARD,FFTW_ESTIMATE);
	memcpy(in,l.data,arrsize);
	fftw_execute(pl);
	memcpy(in,r.data,arrsize);
	fftw_execute(pr);
	//cblas_zsbmv('U',N,0,FFTW_ONE,)
	mult_arr(lout,rout,in,N);
	fftw_execute(pc);
	TSignal<cx_double> s(lout,l.t,N,Signal::SAME);
	fftw_free(in);
	fftw_free(lout);
	fftw_free(rout);
	fftw_destroy_plan(pl);
	fftw_destroy_plan(pr);
	fftw_destroy_plan(pc);
	return s;
}


Signal::TSignal<cx_double> Signal::crosscorr(const Signal::TSignal<cx_double>& l, const Signal::TSignal<cx_double>& r){
	int N = l.len;
	int arrsize = N*sizeof(fftw_complex);

	fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *lout = (fftw_complex*)fftw_malloc(arrsize);
	fftw_complex *rout = (fftw_complex*)fftw_malloc(arrsize);

	fftw_plan pl = fftw_plan_dft_1d(N,in,lout,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan pr = fftw_plan_dft_1d(N,in,rout,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan pc = fftw_plan_dft_1d(N,in,lout,FFTW_BACKWARD,FFTW_ESTIMATE);

	memcpy(in,l.data,arrsize);
	conj_arr(in,N);
	fftw_execute(pl);
	memcpy(in,r.data,arrsize);
	fftw_execute(pr);
  //cblas_zsbmv('U',N,0,FFTW_ONE,...) // can be used for fast element-wise vector multiplication
	mult_arr(lout,rout,in,N);
	fftw_execute(pc);
	TSignal<cx_double> s(lout,l.t,N,Signal::SAME);
	fftw_free(in);
	fftw_free(lout);
	fftw_free(rout);
	fftw_destroy_plan(pl);
	fftw_destroy_plan(pr);
	fftw_destroy_plan(pc);
	return s;
}


Signal::TSignal<cx_double> Signal::autocorr(const Signal::TSignal<cx_double>& f){
	return crosscorr(f, f);
}


Signal::TSignal<double>& Signal::fftshift(Signal::TSignal<double>& f){
	size_t left = f.len/2;
	size_t right = f.len-left;
	size_t n_bytes = sizeof(double);

	double *tmp = (double*)malloc(n_bytes*left);
	memmove(tmp,f.data+right,left*n_bytes);
	memmove(f.data+left,f.data,right*n_bytes);
	memmove(f.data,tmp,left*n_bytes);
	return f;
}
