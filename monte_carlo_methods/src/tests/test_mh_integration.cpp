//created by John Stanco 7.13.18

#include "../include/help.hpp"
//#include "../include/timer.hpp"

typedef double ( *function )( const arma::vec & );



double propose_new( const double mu, const double var )
{
	return 1;
}	


double accept( const double f, const double f0 )
{
	if ( f >= f0 ) return 1;
	if( f / f0 > fRand( 0, 1 ) ) return 1;
	return 0;
}


double harmonic_estimator_integral( const function fn, const arma::vec & xmin, const arma::vec & xmax )
{

	arma::vec dx 	= xmax - xmin;
	arma::vec x0 	= uniDistVec( xmin, xmax );
	//double wid 		= std::abs( xmax - xmin ) / 10;
	double f0 		= fn( x0 );
	double I 		= 0;

	for( size_t i = 0; i < 10000; i++ )
	{
		arma::vec x = uniDistVec( xmin, xmax );
		double f = fn( x );
		if( accept( f, f0 ) ){
			x0 = x;
			f0 = f;
		} 
	}

	//mc loop
	size_t iter = 1e6;
	for( size_t i = 0; i < iter; i++ )
	{
		arma::vec x = uniDistVec( xmin, xmax );
		double f = fn( x );

		if( accept( f, f0 ) ){
			x0 = x;
			f0 = f;
		} 
		
		I += pow( f0, -1 );
	}
	return pow( I / ( arma::prod( dx ) * iter ) , -1 );
}



double basic_mc_integration( const function fn, const arma::vec & xmin, const arma::vec & xmax )
{
	size_t iter = 1e6;
	double I = 0;
	for( size_t i = 0; i < iter; i++ )
	{
		arma::vec x = uniDistVec( xmin, xmax );
		I += fn( x );
	}
	return I * prod( xmax - xmin ) / iter;
}


double importance_sampling_mc( const function fn, const arma::vec & xmin, const arma::vec & xmax, const function p )
{
	arma::vec dx 	= xmax - xmin;
	arma::vec x0 	= uniDistVec( xmin, xmax );
	//double wid 		= std::abs( xmax - xmin ) / 10;
	double p0 		= p( x0 );
	double I 		= 0;

	size_t iter = 1e5;
	size_t burn = 10000;
	size_t lag 	= 2;

	//do for positive and negative part...

	for( size_t i = 0; i < burn; i++ )
	{
		arma::vec x = uniDistVec( xmin, xmax );
		double px = p( x );
		if( accept( px, p0 ) ){
			x0 = x;
			p0 = px;
		} 
	}

	double f = fn( x0 );
	
	for( size_t i = 0; i < iter; i++ )
	{
		arma::vec x = uniDistVec( xmin, xmax );
		double px = p( x );
		if( accept( px, p0 ) ){
			x0 = x;
			p0 = px;
			f  = fn( x0 );
		} 
		if( i % lag == 0 ){
			//x0.print();
			I += f / p0;
		}
	}
	return I / iter * lag;
}



double multi_gaussian( const arma::vec & x )
{
	double std = .5;
	return pow( 2 * pi, -.5 * x.n_elem ) * exp( -.5 * dot( x, x ) / ( std * std ) ) / std / .9539; 
}


double multi_gaussian_2( const arma::vec & x )
{
	double std = .54;
	return pow( 2 * pi, -.5 * x.n_elem ) * exp( -.5 * dot( x, x ) / ( std * std ) ) / std; 
}


double one( const arma::vec & x ) 
{
	return .5;
}


double lin( const arma::vec & x )
{
	return arma::sum( x );
}


double test_convex( const arma::vec & x )
{
	return 1 - pow( sqrt( dot( x, x ) ), 1.6 );
}


int main(){

	arma::vec xmin, xmax;

	xmin << -1 << arma::endr;
	xmax << 1 << arma::endr;

	//std::cout << harmonic_estimator_integral( lin, xmin, xmax ) << std::endl;
	std::cout << basic_mc_integration( multi_gaussian_2, xmin, xmax ) << std::endl;
	std::cout << importance_sampling_mc( multi_gaussian_2, xmin, xmax, multi_gaussian ) << std::endl;
	std::cout << importance_sampling_mc( multi_gaussian_2, xmin, xmax, one ) << std::endl;
	return 1;
}

//store a disconnected correlation matrix...