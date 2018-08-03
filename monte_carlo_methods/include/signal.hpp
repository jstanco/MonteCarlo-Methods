//Created by John Stanco on 2/5/18

#ifndef TSignal_H
#define TSignal_H

#include "timer.hpp"
#include "mc_methods.hpp"


typedef std::complex<double> cx_double;
extern double pi;
//double pi = std::acos(-1);

namespace Signal{

	/// helper functions for use with fftw types

	inline void assign(fftw_complex z1, const fftw_complex z2){
		z1[0]=z2[0];
		z1[1]=z2[1];
	}

	inline void assign_conj(fftw_complex z1, const fftw_complex z2){
		z1[0]=z2[0];
		z1[1]=-z2[1];
	}

	inline void multiply(fftw_complex l, fftw_complex r, fftw_complex rslt){
		rslt[0] = l[0]*r[0]-l[1]*r[1];
		rslt[1] = l[0]*r[1]+l[1]*r[0];
	}

	inline void mult_arr(fftw_complex l[], fftw_complex r[], fftw_complex rslt[], const size_t len){
		for( size_t i=0;i<len;i++){
			multiply(l[i],r[i],rslt[i]);
		}
	}

	inline void conj(fftw_complex z){ z[1]*=-1; }

	inline void conj_arr(fftw_complex f[], const size_t len){
		for( size_t i=0;i<len;i++){ conj(f[i]); }
	}

	enum fft_id { SAME, DUAL };
	
	template<class T>
	class TSignal{
	private:
		T *data;
		double *t;
		size_t len;
		size_t cap;

		friend void swap(TSignal& a, TSignal& b){
			std::swap(a.data, b.data);
			std::swap(a.t, b.t);
			std::swap(a.len, b.len);
			std::swap(a.cap, b.cap);
		}

		void populate_t(const double tmin, const double tmax);
		void populate_data(T(*f)(double));

		TSignal(double const *_t, const size_t n);
		TSignal(fftw_complex const *fftw_data, double const *_t, const size_t n, fft_id id);
		TSignal(double const *fftw_data, double const *_t, const size_t n, fft_id id);
	public:
		TSignal() : len(0), data(0), t(0) {}
		TSignal(const size_t n) : len(n), data(new T[n]), t(new double[n]) { populate_t(0,n); }
		TSignal(const size_t n, const double tmin, const double tmax);
		TSignal(const T _data[], double const *_t, const size_t n);
		TSignal(T(*f)(double), double const *_t, const size_t n);
		TSignal(T(*f)(double), const double tmin, const double tmax, const size_t n);
		TSignal(const TSignal& other);
		~TSignal(){delete [] data; delete [] t;}

		TSignal& operator=(TSignal other){ swap(*this, other); return *this; }
		TSignal operator*(const TSignal& other);

		/// TODO:
		/// add operators that use BLAS for efficient operations
		/// template specialize double and complex TSignals

		T& operator()(const size_t index) { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		T& operator()(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		double time(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return t[index]; }
		size_t n_sample() const { return len; }
		void record(const T& datum, const size_t index){ data[index]=datum; }

		friend T mean(const TSignal&);
		friend T var(const TSignal&);
	};


	template<>
	class TSignal<cx_double>{
	private:
		cx_double *data;
		double *t;
		size_t len;
		size_t cap;

		friend void swap(TSignal& a, TSignal& b){
			std::swap(a.data, b.data);
			std::swap(a.t, b.t);
			std::swap(a.len, b.len);
			std::swap(a.cap, b.cap);
		}
		void populate_t(const double tmin, const double tmax);
		void populate_data(cx_double(*f)(double));

		TSignal(fftw_complex const *fftw_data, double const *_t, const size_t n, fft_id id);
		TSignal(double const *_t, size_t n);

	public:
		TSignal() : len(0), data(0), t(0) {}
		TSignal(const size_t n) : len(n), data(new cx_double[n]), t(new double[n]) { populate_t(0,n); }
		TSignal(const cx_double _data[], double const *_t, const size_t n);
		TSignal(cx_double(*f)(double), double const *_t, const size_t n);
		TSignal(cx_double(*f)(double), const double tmin, const double tmax, const size_t n);
		TSignal(const TSignal& other);
		~TSignal(){delete [] data; delete [] t;}

		TSignal& operator=(TSignal other){ swap(*this, other); return *this; }
		TSignal operator*(const TSignal& other);

		/// TODO:
		/// add operators that use BLAS for efficient operations
		/// template specialize double and complex TSignals

		cx_double& operator()(const size_t index) { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		cx_double& operator()(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		double time(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return t[index]; }
		size_t n_sample() const { return len; }
		void record(const cx_double& datum, const size_t index){ data[index]=datum; }

		friend TSignal<cx_double> convolve(const TSignal<cx_double>& l, const TSignal<cx_double>& r);
		friend TSignal<cx_double> crosscorr(const TSignal<cx_double>& l, const TSignal<cx_double>& r);
		friend TSignal<cx_double> autocorr(const TSignal<cx_double>& f);
		friend TSignal<cx_double>& fftshift(TSignal<cx_double>& f);
		friend cx_double mean(const TSignal<cx_double>& f);
		friend cx_double var(const TSignal<cx_double>& f);
		friend TSignal<cx_double> fft(const TSignal<cx_double>& f);
		friend TSignal<cx_double> ifft(const TSignal<cx_double>& f);
		friend TSignal<cx_double> fft_r2c(const TSignal<double>& f);
		friend TSignal<double> ifft_c2r(const TSignal<cx_double>& f);
	};


	template<>
	class TSignal<double>{
	private:
		double *data;
		double *t;
		size_t len;
		size_t cap;

		friend void swap(TSignal& a, TSignal& b){
			std::swap(a.data, b.data);
			std::swap(a.t, b.t);
			std::swap(a.len, b.len);
			std::swap(a.cap, b.cap);
		}
		void populate_t(const double tmin, const double tmax);
		void populate_data(double(*f)(double));

		TSignal(double const*, double const*_t, const size_t n, fft_id id);
		TSignal(double const *_t, size_t n);

	public:
		TSignal() : len(0), data(0), t(0) {}
		TSignal(const size_t n) : len(n), data(new double[n]), t(new double[n]) { populate_t(0,n); }
		TSignal(const size_t n, const double tmin, const double tmax);
		TSignal(const double _data[], double const *_t, const size_t n);
		TSignal(double(*f)(double), double const *_t, const size_t n);
		TSignal(double(*f)(double), const double tmin, const double tmax, size_t n);
		TSignal(const TSignal& other);
		~TSignal(){ delete [] data; delete [] t; }

		TSignal& operator=(TSignal other){ swap(*this, other); return *this; }
		TSignal operator*(const TSignal& other);

		/// TODO:
		/// add operators that use BLAS for efficient operations
		/// template specialize double and complex TSignals

		double& operator()(const size_t index) { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		double& operator()(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return data[index]; }
		double time(const size_t index) const { if(index>=len){ throw "Index out of bounds."; } return t[index]; }
		size_t n_sample() const { return len; }
		void record(const double& datum, const size_t index){ data[index]=datum; }

		friend TSignal<double> convolve(const TSignal<double>& l, const TSignal<double>& r);
		friend TSignal<double> crosscorr(const TSignal<double>& l, const TSignal<double>& r);
		friend TSignal<double> autocorr(const TSignal<double>& f);
		friend TSignal<double>& fftshift(TSignal<double>& f);
		friend double mean(const TSignal<double>& f);
		friend double var(const TSignal<double>& f);
		friend TSignal<double> fft(const TSignal<double>& f);
		friend TSignal<double> ifft(const TSignal<double>& f);
		friend TSignal<cx_double> fft_r2c(const TSignal<double>& f);
		friend TSignal<double> ifft_c2r(const TSignal<cx_double>& f);
	};


	//////////////////// IMPLEMEMTATION ///////////////////////


	inline
	void TSignal<double>::populate_t(const double tmin, const double tmax){
		double t_curr = tmin;
		double dt = (tmax-tmin)/(len-1);
		for(size_t i=0;i<len;i++){
			t[i] = t_curr;
			t_curr+=dt;
		}
	}


	inline
	void TSignal<cx_double>::populate_t(const double tmin, const double tmax){
		double t_curr = tmin;
		double dt = (tmax-tmin)/(len-1);
		for(size_t i=0;i<len;i++){
			t[i] = t_curr;
			t_curr+=dt;
		}
	}


	//template<>
	inline
	TSignal<cx_double>::TSignal(fftw_complex const *fftw_data, double const *_t, const size_t n, fft_id id) : 
	len(n), data(new cx_double[n]), t(new double[n]){
		memcpy(data, fftw_data, len*sizeof(cx_double));
		if(id == SAME){ memcpy(t, _t, len*sizeof(double)); }
		else {
			double f_n = -pow(_t[len-1]-_t[0], -1)*len*pi; //nyquist frequency
			populate_t(-f_n, f_n);	
		}
	}


	//template<>
	inline
	TSignal<double>::TSignal(double const *fftw_data, double const *_t, const size_t n, fft_id id) : 
	len(n),cap(n),data(new double[n]),t(new double[n]){
		memcpy(data, fftw_data, len*sizeof(double));
		if(id == SAME){ memcpy(t, _t, len*sizeof(double)); }
		else {
			double f_n = -pow(_t[len-1]-_t[0], -1)*len*pi; //nyquist frequency
			populate_t(-f_n, f_n);	
		}
	}


	inline
	TSignal<double>::TSignal(const size_t n, const double tmin, const double tmax) :
	len(n), cap(n), data(new double[n]), t(new double[n])	{
		populate_t(tmin, tmax);
	}


	inline
	TSignal<double>::TSignal(TSignal<double> const& other) :
	len(other.len), cap(other.cap), data(new double[other.len]), t(new double[other.len]){
		memcpy(data, other.data, len*sizeof(double));
		memcpy(t, other.t, len*sizeof(double));
	}


	template<class T>
	void TSignal<T>::populate_t(const double tmin, const double tmax){
		double t_curr = tmin;
		double dt = (tmax-tmin)/(len-1);
		for(size_t i=0;i<len;i++){
			t[i] = t_curr;
			t_curr+=dt;
		}
	}


	template<class T>
	void TSignal<T>::populate_data(T(*f)(double)){
		for(size_t i=0; i<len;i++){ data[i] = f(t[i]); }
	}

	
	/// uninitialized TSignal
	template<class T>
	TSignal<T>::TSignal(double const *_t, const size_t n) : 
	len(n),cap(n),data(new T[n]),t(new double[n]) {
		memcpy(t, _t, len*sizeof(double));
	}


	template<class T>
	TSignal<T>::TSignal(const size_t n, const double tmin, const double tmax) :
	len(n),cap(n),data(new T[n],t(new double[n])) {
		populate_t(tmin,tmax);
	}


	template <class T>
	TSignal<T>::TSignal(T const *_data, double const *_t, const size_t n) :
	len(n), cap(len), data(new T(n)), t(new double[n]) {
		memcpy(data, _data, len*sizeof(T));
		memcpy(t, _t, len*sizeof(double));
	}


	template <class T>
	TSignal<T>::TSignal(T(*f)(double), double const *_t, const size_t n) : 
	len(n), data(new T[n]), t(new double[n]) {
		memcpy(t, _t, len*sizeof(double));
		populate_data(f);
	}


	template <class T>
	TSignal<T>::TSignal(T(*f)(double), const double tmin, const double tmax, const size_t n) : 
	len(n), data(new T[n]), t(new double[n]) {
		populate_t(tmin,tmax);
		populate_data(f);
	}


	template <class T>
	TSignal<T>::TSignal(const TSignal<T>& other) :
	len(other.len), data(new T[other.len]), t(new double[other.len]) {
		memcpy(data, other.data, len*sizeof(T));
		memcpy(t, other.t, len*sizeof(double));
	}


	template<class T>
	TSignal<T> TSignal<T>::operator*(const TSignal& other){
		if(len != other.len){ throw "TSignals must be of same size"; }
		TSignal tmp(t,len);
		for(size_t i=0;i<len;i++){ tmp.data[i]=data[i]*other.data[i]; }
		return tmp;			
	}

	void print(const TSignal<cx_double>&);
	void print_real(const TSignal<cx_double>&, const char*);
	void print_imag(const TSignal<cx_double>&, const char*);
	void print_real(const TSignal<cx_double>&, const std::string&);
	void print_imag(const TSignal<cx_double>&, const std::string&);
	void print(const TSignal<double>&);
	void print(const TSignal<double>&, const char*);
	void print(const TSignal<double>&, const std::string&);
  
	typedef TSignal<cx_double> cx_signal;
	typedef TSignal<double> signal;
}

#endif /* TSignal_H */