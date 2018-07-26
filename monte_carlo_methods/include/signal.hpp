//Created by John Stanco on 2/5/18

#ifndef signal_h
#define signal_h

#include "timer.hpp"
#include "mc_methods.hpp"


typedef std::complex<double> cx_double;
extern double pi;
//double pi = std::acos(-1);

namespace Signal {

	/// helper functions for use with fftw types

	void assign(fftw_complex z1, const fftw_complex z2){
		z1[0]=z2[0];
		z1[1]=z2[1];
	}

	void assign_conj(fftw_complex z1, const fftw_complex z2){
		z1[0]=z2[0];
		z1[1]=-z2[1];
	}

	void multiply(fftw_complex l, fftw_complex r, fftw_complex rslt){
		rslt[0] = l[0]*r[0]-l[1]*r[1];
		rslt[1] = l[0]*r[1]+l[1]*r[0];
	}

	void mult_arr(fftw_complex l[], fftw_complex r[], fftw_complex rslt[], const size_t len){
		for( size_t i=0;i<len;i++){
			multiply(l[i],r[i],rslt[i]);
		}
	}

	void conj(fftw_complex z){ z[1]*=-1; }

	void conj_arr(fftw_complex f[], const size_t len){
		for( size_t i=0;i<len;i++){ conj(f[i]); }
	}

	enum fft_id{
		SAME=1,
		DUAL=2
	};
	
	template<class T>
	class t_signal{
	private:
		T *data;
		double *t;
		size_t len;
		size_t cap;

		friend void swap(t_signal& a, t_signal& b){
			std::swap(a.data, b.data);
			std::swap(a.t, b.t);
			std::swap(a.len, b.len);
			std::swap(a.cap, b.cap);
		}
		void populate_t(double tmin, double tmax);
		void populate_data(T(*f)(double));

		t_signal<cx_double>(fftw_complex fftw_data[], const double _t[], const size_t n, fft_id id);
		t_signal<double>(double const*, double const*_t, const size_t n, fft_id id);
		t_signal(const double _t[], size_t n);

	public:
		t_signal() : len(0), data(0), t(0) {}
		t_signal(size_t n) : len(n), data(new T[n]), t(new double[n]) { populate_t(0,n); }
		t_signal(const T _data[], const double _t[], const size_t n);
		t_signal(T(*f)(double), const double _t[], size_t n);
		t_signal(T(*f)(double), const double tmin, const double tmax, size_t n);
		t_signal(const t_signal& other);
		~t_signal(){delete [] data; delete [] t;}

		t_signal& operator=(t_signal other){ swap(*this, other); return *this; }
		t_signal operator*(const t_signal& other);

		/// TODO:
		/// add operators that use BLAS for efficient operations
		/// template specialize double and complex signals

		T operator()(const size_t index) const { if(index>=len){ throw "Index our of bounds."; } return data[index]; }
		double time(const size_t index) const { if(index>=len){ throw "Index our of bounds."; } return t[index]; }
		size_t n_sample() const { return len; }
		void record(const T& datum, const size_t index){ data[index]=datum; }

		friend t_signal<cx_double> fft(const t_signal<cx_double>& f);
		friend t_signal<cx_double> ifft(const t_signal<cx_double>& f);
		friend t_signal<double> fft_r2r(const t_signal<double>& f);
		friend t_signal<double> ifft_r2r(const t_signal<double>& f);
		friend t_signal<cx_double> fft_r2c(const t_signal<double>& f);
		friend t_signal<double> ifft_c2r(const t_signal<cx_double>& f);
		friend t_signal convolve(const t_signal& l, const t_signal& r);
		friend t_signal crosscorr(const t_signal& l, const t_signal& r);
		friend t_signal autocorr(const t_signal& f);
		friend t_signal autocorr_slow(const t_signal& f);
		friend t_signal& fftshift(t_signal& f);
		friend T mean(const t_signal& f);
		friend T var(const t_signal& f);
	};




	//////////////////// IMPLEMEMTATION ///////////////////////


	template<class T>
	void t_signal<T>::populate_t(double tmin, double tmax){
		double t_curr = tmin;
		double dt = (tmax-tmin)/(len-1);
		for(size_t i=0; i<len;i++){
			t[i] = t_curr;
			t_curr+=dt;
		}
	}


	template<class T>
	void t_signal<T>::populate_data(T(*f)(double)){
		for(size_t i=0; i<len;i++){ data[i] = f(t[i]); }
	}

	
	/// uninitialized signal
	template<class T>
	t_signal<T>::t_signal(const double _t[], size_t n) : 
	len(n),cap(n),data(new T[n]),t(new double[n]) {
		memcpy(t, _t, len*sizeof(double));
	}


	template<>
	t_signal<cx_double>::t_signal(fftw_complex fftw_data[], const double _t[], size_t n, fft_id id) : 
	len(n), data(new cx_double[n]), t(new double[n]){
		memcpy(data, fftw_data, len*sizeof(cx_double));
		if(id == SAME){ memcpy(t, _t, len*sizeof(double)); }
			else {double f_n = -pow(_t[len-1]-_t[0], -1)*len*pi; //nyquist frequency
								populate_t(-f_n, f_n);	
		}
	}


	template <class T>
	t_signal<T>::t_signal(const T _data[], const double _t[], const size_t n) :
	len(n), cap(len), data(new T(n)), t(new double[n]) {
		memcpy(data, _data, len*sizeof(T));
		memcpy(t, _t, len*sizeof(double));
	}


	template <class T>
	t_signal<T>::t_signal(T(*f)(double), const double _t[], size_t n) : 
	len(n), data(new T[n]), t(new double[n]) {
		memcpy(t, _t, len*sizeof(double));
		populate_data(f);
	}


	template <class T>
	t_signal<T>::t_signal(T(*f)(double), const double tmin, const double tmax, size_t n) : 
	len(n), data(new T[n]), t(new double[n]) {
		populate_t(tmin,tmax);
		populate_data(f);
	}


	template <class T>
	t_signal<T>::t_signal(const t_signal<T>& other) :
	len(other.len), data(new T[other.len]), t(new double[other.len]) {
		memcpy(data, other.data, len*sizeof(T));
		memcpy(t, other.t, len*sizeof(double));
	}


	template<class T>
	t_signal<T> t_signal<T>::operator*(const t_signal& other){
		if(len != other.len){ throw "Signals must be of same size"; }
		t_signal tmp(t,len);
		for(size_t i=0;i<len;i++){ tmp.data[i]=data[i]*other.data[i]; }
		return tmp;			
	}


	/////////////////////// REAL SIGNAL ////////////////////////


	double mean(const t_signal<double>& f){
		double mu = 0;
		for(size_t i=0;i<f.len;i++){
			mu += f(i);
		}
		return mu / f.len;
	}


	double var(const t_signal<double>& f){
		double sigma = 0;
		double mu = mean(f);
		for(size_t i=0;i<f.len;i++){
			sigma += pow(f(i) - mu,2);
		}
		return sigma / f.len;
	}


	template<>
	t_signal<double>::t_signal(double const *fftw_data, double const *_t, const size_t n, fft_id id) : 
	len(n),cap(n),data(new double[n]),t(new double[n]){
		memcpy(data, fftw_data, len*sizeof(double));
		if(id == SAME){ memcpy(t, _t, len*sizeof(double)); }
		else {
				double f_n = -pow(_t[len-1]-_t[0], -1)*len*pi; //nyquist frequency
				populate_t(-f_n, f_n);	
		}
	}


	void print(const t_signal<double>& f){
		for(size_t i=0;i<f.n_sample();i++){
			printf("t: %lf\tf: %lf\n", f.time(i),f(i));
		}
	}


	void print(const t_signal<double> &f, const char* filename){
		FILE *pFile = fopen(filename, "w");
		if(pFile){ for(size_t i=0;i<f.n_sample();i++){
				fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i)); }
		}
	}

	void print(const t_signal<double> &f, const std::string& filename){
		FILE *pFile = fopen(filename.c_str(), "w");
		if(pFile){ for(size_t i=0;i<f.n_sample();i++){
				fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i)); }
		}
	}


	/// DCT-II - input data must be real and symmetric
	t_signal<double> fft_r2r(const t_signal<double>& f){
		int arrsize = f.len*sizeof(double);
		double *in = (double*)fftw_malloc(arrsize);
		double *out = (double*)fftw_malloc(arrsize);
		fftw_plan p = fftw_plan_r2r_1d(f.len, in, out, FFTW_REDFT10,FFTW_ESTIMATE);
		memcpy(in,f.data,arrsize);
		fftw_execute(p);
		t_signal<double> s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}


	t_signal<double> ifft_r2r(const t_signal<double>& f){
		int arrsize = f.len*sizeof(double);
		double *in = (double*)fftw_malloc(arrsize);
		double *out = (double*)fftw_malloc(arrsize);	//exploits fact that fft of real function is hermitian
		fftw_plan p = fftw_plan_r2r_1d(f.len, in, out,FFTW_REDFT01,FFTW_ESTIMATE);
		memcpy(in,f.data,arrsize);
		fftw_execute(p);
		t_signal<double> s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}


	/// IDCT-II - input data must be real and symmetric
	t_signal<cx_double> fft_r2c(const t_signal<double>& f){
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

		t_signal<cx_double> s(rslt,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		fftw_free(rslt);
		return s;
	}


	t_signal<double> ifft_c2r(const t_signal<cx_double>& f){
		int in_arrsize = (f.len/2+1)*sizeof(fftw_complex);
		int out_arrsize = f.len*sizeof(double);
		fftw_complex *in = (fftw_complex*)fftw_malloc(in_arrsize);
		double *out = (double*)fftw_malloc(out_arrsize);
		fftw_plan p = fftw_plan_dft_c2r_1d(f.len, in, out, FFTW_ESTIMATE);
		memcpy(in,f.data,in_arrsize);
		fftw_execute(p);
		t_signal<double> s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}


	t_signal<double> convolve(const t_signal<double>& l, const t_signal<double>& r){
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
		t_signal<double> s(in,l.t,n,SAME);
		fftw_free(in);
		fftw_free(lout);
		fftw_free(rout);
		fftw_free(rslt);
		fftw_destroy_plan(pl);
		fftw_destroy_plan(pr);
		fftw_destroy_plan(pc);
		return s;
	}


	t_signal<double> crosscorr(const t_signal<double>& l, const t_signal<double>& r){
		return convolve(l,r);
	}

	t_signal<double> autocorr(const t_signal<double>& f){
		return crosscorr(f,f);
	}

	///////////////////// COMPLEX SIGNAL ///////////////////////



	void print(const t_signal<cx_double>& f){
		for(size_t i=0;i<f.n_sample();i++){
			printf("t: %lf\tre: %lf\tim: %lf\n", f.time(i),f(i).real(),f(i).imag());
		}
	}

	void print_real(const t_signal<cx_double> &f, const char* filename){
		FILE *pFile = fopen(filename, "w");
		if(pFile){ for(size_t i=0;i<f.n_sample();i++){
				fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i).real()); }
		}
	}

	void print_imag(const t_signal<cx_double> &f, const char* filename){
		FILE *pFile = fopen(filename, "w");
		if(pFile){ for(size_t i=0;i<f.n_sample();i++){
				fprintf(pFile,"%lf\t\t%lf\n", f.time(i), f(i).imag()); }
		}
	}


	t_signal<cx_double> fft(const t_signal<cx_double>& f){
		int arrsize = sizeof(fftw_complex)*f.len;
		fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
		fftw_complex *out = (fftw_complex*)fftw_malloc(arrsize);
		fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
		memcpy(in,f.data,arrsize);
		fftw_execute(p);
		t_signal<cx_double> s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}


	t_signal<cx_double> ifft(const t_signal<cx_double>& f){
		int arrsize = sizeof(fftw_complex)*f.len;
		fftw_complex *in = (fftw_complex*)fftw_malloc(arrsize);
		fftw_complex *out = (fftw_complex*)fftw_malloc(arrsize);
		fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
		memcpy(in,f.data,arrsize);
		fftw_execute(p);
		t_signal<cx_double> s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}


	/// assume same length for now
	t_signal<cx_double> convolve(const t_signal<cx_double>& l, const t_signal<cx_double>& r){
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
		t_signal<cx_double> s(lout,l.t,N,SAME);
		fftw_free(in);
		fftw_free(lout);
		fftw_free(rout);
		fftw_destroy_plan(pl);
		fftw_destroy_plan(pr);
		fftw_destroy_plan(pc);
		return s;
	}


	t_signal<cx_double> crosscorr(const t_signal<cx_double>& l, const t_signal<cx_double>& r){
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
		t_signal<cx_double> s(lout,l.t,N,SAME);
		fftw_free(in);
		fftw_free(lout);
		fftw_free(rout);
		fftw_destroy_plan(pl);
		fftw_destroy_plan(pr);
		fftw_destroy_plan(pc);
		return s;
	}


	t_signal<cx_double> autocorr(const t_signal<cx_double>& f){
		return crosscorr(f, f);
	}


	t_signal<double>& fftshift(t_signal<double>& f){
		size_t left = f.len/2;
		size_t right = f.len-left;
		size_t n_bytes = sizeof(double);

		double *tmp = (double*)malloc(n_bytes*left);
		memmove(tmp,f.data+right,left*n_bytes);
		memmove(f.data+left,f.data,right*n_bytes);
		memmove(f.data,tmp,left*n_bytes);
		return f;
	}

	t_signal<double> autocorr_slow(const t_signal<double>& f){
		t_signal<double> s(f.t, f.len);
		for(size_t i=0;i<f.len;i++){
			s.record(0,i);
			for(size_t j=0;j<f.len;j++){
				s.data[i]+=f(j)*f((i+j)%f.len);
			}
		}
		return s;
	}
  
	typedef t_signal<cx_double> cx_signal;
	typedef t_signal<double> signal;
}

#endif /* signal_h */