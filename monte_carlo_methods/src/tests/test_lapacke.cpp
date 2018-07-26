//created by John Stanco 7.19.18



#include "../include/help.hpp"
#include "../include/timer.hpp"

#define SPACE printf("\n")

/*
extern "C" void dgemm_(char *transa, char *transb,
                       int *m, int *n, int *k, double *alpha, double *a, int *lda,
                       double *b, int *ldb, double *beta, double *c, int *ldc );
*/

void rand_mat_colmajor(size_t m, size_t n, double *a, size_t lda) {
	size_t i,j;
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			a[i+lda*j] = fRand(0, 1);
		}
	}
}


void print_mat(size_t m, size_t n, double *a, size_t lda) {
	size_t i,j;
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			printf("%lf ",a[i+lda*j]);
		}
		printf("\n");
	}
}


void zeros(double *arr, const size_t len){
	for(size_t i=0;i<len;i++){arr[i]=0;}
}


template <class T>
arma::Mat<T> arr_to_arma(size_t m, size_t n, T *a, size_t lda){
	arma::Mat<T> A(m, n);
	size_t i,j;
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			A(i,j) = a[i+lda*j];
		}
	}
	return A;
}


void test_arr_to_arma(){
	int m,n,lda;
	double a[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};

	m = 5;
	n = 3;
	lda = 5;

	arma::mat A = arr_to_arma(m,n,a,lda);
	print_mat(m,n,a,lda);
	A.print();
}


class test_dmat{
private:
	double *data;
	size_t n_row;
	size_t n_col;
	char 	*trans;

	void zeros(double* c, const size_t len){ for(size_t i=0;i<len;i++){ c[i] = 0; } }
	void copyarr(double* a, const double*b, const size_t len){ for(size_t i=0;i<len;i++){ a[i] = b[i]; } }
public:
	test_dmat(double a[], size_t m, size_t n) : data(new double[m*n]),n_row(m),n_col(n),trans("n") { copyarr(data,a,m*n); }
	test_dmat(const test_dmat& other) : data(new double[other.n_row*other.n_col]), 
	n_row(other.n_row),n_col(other.n_col),trans(other.trans) { copyarr(data,other.data,n_row*n_col); }
	~test_dmat(){ delete[] data; }

	test_dmat operator*(test_dmat& other){
		int m,n,k,lda,ldb,ldc,clen;
		double alpha, beta;
		alpha = 1.;
		beta = 1.;

		if(n_col != other.n_row) {throw "Matrix dimensions mismatched."; }

		m = n_row;
		n = other.n_col;
		k = n_col;
		lda = m;
		ldb = k;
		ldc = m;
		clen = m*n;
		alpha = 1.;
		beta = 1.;

		double c[clen];
		zeros(c, clen);

		dgemm_(trans,other.trans,&m,&n,&k,&alpha,data,&lda,other.data,&ldb,&beta,c,&ldc);
		return test_dmat(c,m,n);
	}
	void print(){ print_mat(n_row,n_col,data,n_row); }
};


void test_dgemm() {
	int info,m,n,k,lda,ldb,ldc;
	double one,alpha,beta;

	double a[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};
	double b[3*2] = {2,1,0,1,1,2};
	double c[5*2] = {0,0,0,0,0,0,0,0,0,0};
  
	m = 5;
	n = 2;
	k = 3;
	lda = m;
	ldb = k;
	ldc = m;
	one = 1.;
	alpha = 2.;
	beta = 3.;

	char transa[] = "n";
	char transb[] = "n";

	timer t;
	t.start();
	dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
	t.stop();
	t.print();

	print_mat(m,k,a,lda);
	printf("\n");
	print_mat(k,n,b,ldb);
	printf("\n");
	print_mat(m,n,c,ldc);

}


void test_dgemm_vs_arma() {
	int info,m,n,k,lda,ldb,ldc;
	double one;

	double a[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};
	double b[3*2] = {2,1,0,1,1,2};
	double c[5*2] = {0,0,0,0,0,0,0,0,0,0};
  double c_[5*2] = {-10,12,14,16,18,-3,14,12,16,16};
  
	m = 5;
	n = 2;
	k = 3;
	lda = m;
	ldb = k;
	ldc = m;
	one = 1.;
	char transa[] = "n";
	char transb[] = "n";
	timer t;

	arma::mat A = arr_to_arma(m,k,a,lda);
	arma::mat B = arr_to_arma(k,n,b,ldb);

	t.start();
	arma::mat C = A * B;
	t.stop();
	t.print();

	t.restart();
	dgemm_(transa,transb,&m,&n,&k,&one,a,&lda,b,&ldb,&one,c,&ldc);
	t.stop();
	t.print();

	C.print();
	print_mat(m,n,c,ldc);
}


void test_dgemm_vs_arma_large() {
	int m,n,k,lda,ldb,ldc;
	double one;

	m = 100;
	k = 100;
	n = 100;

	double a[m*k];
	double b[k*n];
	double c[m*n];

	zeros(c,m*n);

	rand_mat_colmajor(m,k,a,m);
	rand_mat_colmajor(k,n,b,k);

	arma::mat A_arma = arr_to_arma(m,k,a,m);
	arma::mat B_arma = arr_to_arma(k,n,b,k);
  
	lda = m;
	ldb = k;
	ldc = m;
	one = 1.;
	char transa[] = "n";
	char transb[] = "n";
	timer t;

	arma::mat A = arr_to_arma(m,k,a,lda);
	arma::mat B = arr_to_arma(k,n,b,ldb);

	t.start();
	arma::mat C = A * B;
	t.stop();
	t.print();

	t.restart();
	dgemm_(transa,transb,&m,&n,&k,&one,a,&lda,b,&ldb,&one,c,&ldc);
	t.stop();
	t.print();

	//C.print();
	//print_mat(m,n,c,ldc);
}


void test_dmat_class() {
	double a[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};
	double b[3*2] = {2,1,0,1,1,2};

  int m,n,k;
	m = 5;
	n = 2;
	k = 3;

	test_dmat A(a,m,k);
	test_dmat B(b,k,n);

	A.print();
	printf("\n");
	B.print();
	printf("\n");

	timer t;
	t.start();
	test_dmat C(A * B);
	t.stop();
	t.print();
	printf("\n");
	C.print();
}


void test_dmat_vs_arma(){
	double a[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};
	double b[3*2] = {2,1,0,1,1,2};

  int m,n,k;
	m = 5;
	n = 2;
	k = 3;

	test_dmat A(a,m,k);
	test_dmat B(b,k,n);

	arma::mat A_arma = arr_to_arma(m,k,a,m);
	arma::mat B_arma = arr_to_arma(k,n,b,k);

	timer t;

	t.start();
	test_dmat C = A * B;
	t.stop();
	t.print();

	t.restart();
	arma::mat C_arma = A_arma * B_arma;
	t.stop();
	t.print();
}


void test_dmat_vs_arma_large() {
	int m,n,k;

	m = 1000;
	k = 1000;
	n = 1000;

	double *a = new double[m*k];
	double *b = new double[k*n];

	rand_mat_colmajor(m,k,a,m);
	rand_mat_colmajor(k,n,b,k);


	test_dmat A(a,m,k);
	test_dmat B(b,k,n);

	arma::mat A_arma = arr_to_arma(m,k,a,m);
	arma::mat B_arma = arr_to_arma(k,n,b,k);

	timer t;

	t.start();
	test_dmat C = A * B;
	t.stop();
	t.print();

	t.restart();
	arma::mat C_arma = A_arma * B_arma;
	t.stop();
	t.print();

	delete [] a;
	delete [] b;
}


int main() {
	test_dgemm_vs_arma_large();
	return 1;
}