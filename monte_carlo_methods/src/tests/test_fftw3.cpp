//created by John Stanco 7.20.18

#include "../../include/help.hpp"
#include "../../include/timer.hpp"
#define I cx_double(0,1)

#ifndef FFTW_ONE
#define FFTW_ONE {1,0}
#endif /* fftw_ONE */

#ifndef FFTW_ZERO
#define FFTW_ZERO {0,0}
#endif /* fftw_ZERO */

#ifndef CX_ONE
#define CX_ONE cx_double(1,0)
#endif /* cx_ONE */

#ifndef CX_ZERO
#define CX_ZERO cx_double(0,0)
#endif /* cx_ZERO */


//should turn this into a class
//signal...
//this is in particular the manifestation of the complex signal...

void print_fftw_complex(const fftw_complex z){
	printf("real: %lf\timag: %lf\n",z[0],z[1]);
}


void print_fftw_complex(const fftw_complex arr[], const size_t len){
	for(size_t i = 0; i < len; i++){
		printf("r: %lf\ti: %f\n",arr[i][0],arr[i][1]);
	}
}


void print_fftw_real(const fftw_complex arr[], const size_t len, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){
		for(size_t i = 0; i < len; i++){
			fprintf(pFile, "%lu\t%f\n",i,arr[i][0]);
		}
	}
	fclose(pFile);
}


void print_fftw_imag(fftw_complex arr[], size_t len, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){
		for(size_t i = 0; i < len; i++){
			fprintf(pFile, "%lu\t%f\n",i,arr[i][1]);
		}
	}
	fclose(pFile);
}


void print_fftw_real(fftw_complex arr[], double t[], size_t len, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){
		for(size_t i = 0; i < len; i++){
			fprintf(pFile, "%lf\t%f\n",t[i],arr[i][0]);
		}
	}
	fclose(pFile);
}

bool equals(const fftw_complex l, const fftw_complex r){
	return (l[0]==r[0] && l[1]==r[1]);
}


void print_fftw_imag(fftw_complex arr[], double t[], size_t len, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){
		for(size_t i = 0; i < len; i++){
			fprintf(pFile, "%lf\t%f\n",t[i],arr[i][1]);
		}
	}
	fclose(pFile);
}


void multiply(fftw_complex left, fftw_complex right, fftw_complex *rslt){
	*rslt[0] = left[0]*right[0] + left[1]*right[1];
	*rslt[1] = left[0]*right[1] + left[1]*right[0];
}


void multiply(fftw_complex left[], fftw_complex right[], fftw_complex rslt[], size_t len){
	for(size_t i=0;i<len;i++){
		multiply(left[i],right[i],&rslt[i]);
	}
}


void to_fftw_complex(fftw_complex *zfft, const cx_double *zcpp){
	memcpy(zfft,zcpp,sizeof(fftw_complex));
}


void to_fftw_complex(fftw_complex *zfft, const double _Complex *zc){
	memcpy(zfft,zc,sizeof(fftw_complex));
}


void exp(fftw_complex arr[], size_t len, double freq, double tmin, double tmax){
	double t = tmin;
	double w = 2*pi*freq;
	double dt = (tmax-tmin)/(len-1);
	for(size_t i=0;i<len;i++){
		cx_double zcpp = exp(I*w*t);
		to_fftw_complex(&arr[i],&zcpp);
		t+=dt;
	}
}


void gaussian(fftw_complex arr[],size_t len, double wid, double tmin, double tmax){
	double t = tmin;
	double dt = (tmax-tmin)/(len-1);
	for(size_t i=0;i<len;i++){
		cx_double zcpp = exp(-.5*t*t/wid);
		to_fftw_complex(&arr[i],&zcpp);
		t+=dt;
	}
}


void linspace(double ts[], size_t len, double tmin, double tmax){
	double t = tmin;
	double dt = (tmax-tmin)/(len-1);
	for(size_t i=0;i<len;i++){
		ts[i] = t;
		t+=dt;
	}
}


void f_to_signal(cx_double (*f)(double), fftw_complex arr[], double t[], size_t len){
	for(size_t i=0;i<len;i++){
		cx_double z = f(t[i]);
		to_fftw_complex(&arr[i],&z);
	}
}


void test_fft_complex(){
	cx_double zcpp(3,1);
	//double _Complex z1 = 3.0+1.0I;
	fftw_complex zfft;
	memcpy(&zfft,&zcpp,sizeof(fftw_complex));
	print_fftw_complex(zfft);

	double _Complex zc = 3.0+1.0I;
	memcpy(&zfft,&zc,sizeof(fftw_complex));
	print_fftw_complex(zfft);
}


cx_double test_f(double t){
	return exp(I * t);
}

cx_double gaussian(double t){
	return exp(-.5*t*t);
}


cx_double box(double t){
	if(t<-1 || t> 1){ return cx_double(0,0); }
	return cx_double (1,0);
}

cx_double abs_exp(const double t){
	return exp(-std::abs(t));
}


void divide(fftw_complex arr[], const size_t len, const double x){
	for(size_t i = 0; i < len; i++){
		arr[i][0] /= x;
		arr[i][1] /= x;
	}
}


template<class T>
void fftshift(T arr[], const size_t len){
	//shifts 0th frequency to center
	size_t left = len/2;
	size_t right = len - left;
	size_t n_bytes = sizeof(T);
	fftw_complex tmp[left];
	memmove(tmp,arr+right,left*n_bytes);
	memmove(arr+left,arr,right*n_bytes);
	memmove(arr,tmp,left*n_bytes);
	free(tmp);
}


void test_fft_complex_1d(){
	int N = pow(2,6);
	fftw_complex in[N], out[N];
	fftw_plan p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	gaussian(in,N,1,-3,3);
	fftw_execute(p);
	fftw_destroy_plan(p);
	print_fftw_real(out,N,"../data/fftw_1d_real.dat");
	print_fftw_imag(out,N,"../data/fftw_1d_imag.dat");
	fftw_plan p_rev = fftw_plan_dft_1d(N,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p_rev);
	fftw_destroy_plan(p_rev);
	divide(in,N,N);
	print_fftw_real(in,N,"../data/ifftw_1d_real.dat");
	print_fftw_imag(in,N,"../data/ifftw_1d_imag.dat");
	fftw_cleanup();
}


void abs(fftw_complex arr[], size_t len){
	for(size_t i = 0;i<len;i++){
		arr[i][0] = std::abs(arr[i][0]);
		arr[i][1] = std::abs(arr[i][1]);
	}
}


double norm(const fftw_complex z){
	return sqrt(z[0]*z[0] + z[1]*z[1]);
}


void norm(fftw_complex arr[], size_t len){
	for(size_t i = 0;i<len;i++){
		arr[i][0]=norm(arr[i]);
		arr[i][1]=0;
	}
}


class t_signal{
private:
	cx_double *data;
	double *t;
	size_t len;
	size_t cap;

	friend void swap(t_signal& a, t_signal& b){
		std::swap(a.data, b.data);
		std::swap(a.t, b.t);
		std::swap(a.len, b.len);
	}

	void populate_t(double tmin, double tmax){
		double t_curr = tmin;
		double dt = (tmax-tmin)/(len-1);
		for(size_t i=0; i<len;i++){
			t[i] = t_curr;
			t_curr+=dt;
		}
	}

	void populate_data(cx_double(*f)(double)){
		for(size_t i=0; i<len;i++){
			data[i] = f(t[i]);
		}
	}

	enum fft_id{
		SAME,
		DUAL
	};

	t_signal(fftw_complex fftw_data[], const double _t[], size_t n, fft_id id) : 
	len(n), data(new cx_double[n]), t(new double[n]){
		memcpy(data, fftw_data, len*sizeof(cx_double));
			if(id == SAME){ memcpy(t, _t, len*sizeof(double)) }
			else {double f_n = -pow(_t[len-1]-_t[0], -1)*len*pi; //nyquist frequency
								populate_t(-f_n, f_n);	
		}
	}

public:
	t_signal() : len(0), data(0), t(0) {}
	t_signal(size_t n) : len(n), data(new cx_double[n]), t(new double[n]) {}
	t_signal(const cx_double _data[], const double _t[], const size_t n) :
	len(n), cap(len), data(new cx_double(n)), t(new double[n]) {
		memcpy(data, _data, len*sizeof(cx_double));
		memcpy(t, _t, len*sizeof(double));
	}

	t_signal(cx_double(*f)(double), const double _t[], size_t n) : 
	len(n), data(new cx_double[n]), t(new double[n]) {
		memcpy(t, _t, len*sizeof(double));
		populate_data(f);
	}

	t_signal(cx_double(*f)(double), const double tmin, const double tmax, size_t n) : 
	len(n), data(new cx_double[n]), t(new double[n]) {
		populate_t(tmin,tmax);
		populate_data(f);
	}

	t_signal(const t_signal& other) :
	len(other.len), data(new cx_double[other.len]), t(new double[other.len]) {
		memcpy(data, other.data, len*sizeof(cx_double));
		memcpy(t, other.t, len*sizeof(double));
	}

	~t_signal(){delete [] data; delete [] t;}
	t_signal& operator=(t_signal other){
		swap(*this, other);
		return *this;
	}

	void print(){
		for(size_t i=0;i<len;i++){
			printf("t: %lf\tre: %lf\tim: %lf\n", t[i],data[i].real(),data[i].imag());
		}
	}

	friend t_signal fft(const t_signal& f){
		fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*f.len);
		fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*f.len);
		fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
		memcpy(in,f.data,f.len*sizeof(fftw_complex));
		fftw_execute(p);
		t_signal s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}

	friend t_signal ifft(const t_signal& f){
		fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*f.len);
		fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*f.len);
		fftw_plan p = fftw_plan_dft_1d(f.len,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
		memcpy(in,f.data,f.len*sizeof(fftw_complex));
		fftw_execute(p);
		t_signal s(out,f.t,f.len,DUAL);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);
		return s;
	}

	/// assume same length for now
	friend t_signal convolve(const t_signal& l, const t_signal& r){
		/*
		size_t N = l.len;
		if(l.len > r.len){
			//pad left
		} else if (l.len < r.len){
			N = r.len;
			//pad right
		}
	*/
		int N,incx,incy;
		N=l.len;
		incx=1;
		incy=1;

		fftw_complex *in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
		fftw_complex *lout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
		fftw_complex *rout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

		fftw_plan pl = fftw_plan_dft_1d(N,in,lout,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_plan pr = fftw_plan_dft_1d(N,in,rout,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_plan pc = fftw_plan_dft_1d(N,in,lout,FFTW_BACKWARD,FFTW_ESTIMATE);

		memcpy(in,l.data,N*sizeof(fftw_complex));
		fftw_execute(pl);
		memcpy(in,r.data,N*sizeof(fftw_complex));
		fftw_execute(pr);
		//cblas_zsbmv('U',N,0,FFTW_ONE,)
		multiply(lout,rout,in,N);
		fftw_execute(pc);

		t_signal s(lout,l.t,N,SAME);

		fftw_free(in);
		fftw_free(lout);
		fftw_free(rout);
		fftw_destroy_plan(pl);
		fftw_destroy_plan(pr);
		fftw_destroy_plan(pc);
		return s;
	}
};


void test_zsbmv(){

}



bool equals(const fftw_complex l[], const fftw_complex r[], size_t len){
	for(size_t i=0;i<len;i++){
		if(!equals(l[i],r[i])){ return false; }
	}
	return true;
}	


void test_fft_shift(){
	int N = pow(2,6);
	fftw_complex in[N], out[N];
	double t[N], f[N];

	fftw_plan p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	double tmin,tmax,T,fmin,fmax;
	tmin = -4;
	tmax = 4;
	T = tmax-tmin;
	fmin = -pow(T,-1)*N*pi;
	fmax = pow(T,-1)*N*pi;


	linspace(t,N,tmin,tmax);
	linspace(f,N,fmin,fmax);

	f_to_signal(&box,in,t,N);

	fftw_execute(p);
	fftw_destroy_plan(p);

	print_fftw_real(out,f,N,"../data/fftw_1d_real.dat");
	print_fftw_imag(out,f,N,"../data/fftw_1d_imag.dat");

	norm(out,N);
	fftshift(out,N);
	
	print_fftw_real(out,f,N,"../data/fftw_1d_shifted_real.dat");
	print_fftw_imag(out,f,N,"../data/fftw_1d_shifted_imag.dat");

	fftw_cleanup();
}


void conj(fftw_complex arr[], size_t len){
	for(size_t i=0;i<len; i++){arr[i][1]*=-1;}
}


void test_convolution(){
	int N = pow(2, 8);
	fftw_complex l_in[N],l_out[N],r_in[N],r_out[N],rev_in[N],rev_out[N];
	double t[N];

	linspace(t,N,0,10);
	fftw_plan p_l = fftw_plan_dft_1d(N,l_in,l_out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan p_r = fftw_plan_dft_1d(N,r_in,r_out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan p_rev = fftw_plan_dft_1d(N,rev_in,rev_out,FFTW_BACKWARD,FFTW_ESTIMATE);

	f_to_signal(&gaussian,l_in,t,N);
	f_to_signal(&box,r_in,t,N);

	print_fftw_real(l_in,t, N,"../data/fftw_convolution.dat");

	fftw_execute(p_l);
	fftw_execute(p_r);

	multiply(l_out,r_out,rev_in,N);
	fftshift(rev_in,N);
	
	fftw_execute(p_rev);
	divide(rev_out,N,N);
	norm(rev_out,N);
	fftw_destroy_plan(p_l);
	fftw_destroy_plan(p_r);
	fftw_destroy_plan(p_rev);
	fftw_cleanup();
}


void test_signal(){
	t_signal sig(abs_exp,0,2,32);
	sig.print();
	t_signal fft_sig = fft(sig);
	fft_sig.print();
}


void test_time(){
	timer t;
	int N = pow(2,16);
	fftw_complex in[N], out[N];
	fftw_plan p = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	gaussian(in,N,1,-3,3);
	t.start();
	fftw_execute(p);
	t.stop(); t.print();
	fftw_destroy_plan(p);

	fftw_complex *in_alloc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *out_alloc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_plan p_alloc = fftw_plan_dft_1d(N,in_alloc,out_alloc,FFTW_FORWARD,FFTW_ESTIMATE);
	gaussian(in_alloc,N,1,-3,3);
	t.restart();
	fftw_execute(p_alloc);
	t.stop();
	t.print();
	fftw_destroy_plan(p_alloc);
	fftw_free(in_alloc);
	fftw_free(out_alloc);
	
	fftw_cleanup();
}


bool test_equals(){
	int N = pow(2, 4);
	fftw_complex *z1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *z2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	double *t = (double*)malloc(sizeof(double)*N);


	linspace(t,N,-3,3);

	f_to_signal(&gaussian, z1, t, N);
	f_to_signal(&gaussian, z2, t, N);

	print_fftw_complex(z1,N);
	printf("\n");
	print_fftw_complex(z2,N);

	bool r = equals(z1, z2, N);

	free(z1);
	free(z2);
	free(t);
	return r;
}


void test_class_convolution(){
	int N = pow(2, 4);
	fftw_complex l_in[N],l_out[N],r_in[N],r_out[N],rev_in[N],rev_out[N];
	double t[N];

	linspace(t,N,0,10);
	fftw_plan p_l = fftw_plan_dft_1d(N,l_in,l_out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan p_r = fftw_plan_dft_1d(N,r_in,r_out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan p_rev = fftw_plan_dft_1d(N,rev_in,rev_out,FFTW_BACKWARD,FFTW_ESTIMATE);

	f_to_signal(&gaussian,l_in,t,N);
	f_to_signal(&box,r_in,t,N);

	//print_fftw_real(l_in,t, N,"../data/fftw_convolution.dat");

	fftw_execute(p_l);
	fftw_execute(p_r);

	multiply(l_out,r_out,rev_in,N);
	//fftshift(rev_in,N);
	fftw_execute(p_rev);
	//divide(rev_out,N,N);
	//norm(rev_out,N);


	t_signal f(&gaussian, t, N);
	t_signal g(&box,t,N);
	t_signal fg = convolve(f,g);

	fg.print();
	print_fftw_complex(rev_out,N);

	fftw_destroy_plan(p_l);
	fftw_destroy_plan(p_r);
	fftw_destroy_plan(p_rev);
	fftw_cleanup();
}


void test_time_of_class(){
	timer t;
	int N = pow(2, 16);

	t.start();
	fftw_complex *in_alloc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *out_alloc = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
	fftw_plan p_alloc = fftw_plan_dft_1d(N,in_alloc,out_alloc,FFTW_FORWARD,FFTW_ESTIMATE);

	double ts[N];
	linspace(ts,N,-3,3);
	f_to_signal(&gaussian,in_alloc,ts,N);

	
	fftw_execute(p_alloc);
	t.stop();
	t.print();

	t_signal sig(gaussian,-3,3,N);
	t.restart();
	t_signal sig_fft = fft(sig);
	t.stop();
	t.print();

	fftw_destroy_plan(p_alloc);

	fftw_free(in_alloc);
	fftw_free(out_alloc);
	fftw_cleanup();
}



int main() {
	test_class_convolution();
	//std::cout << test_equals() << std::endl;
	//test_time_of_class();
	//test_time();
	//test_signal();
	//test_convolution();
	//test_fft_shift();
	//test_fft_complex_1d();
}