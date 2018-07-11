//created by John Stanco 7.10.18


#include "../../Signal/include/signal"


//This will provide post-processing for the Monte-Carlo...

// The point of Monte-Carlo is to create a distribution...
// The monte-carlo samples from a distribution...
// Actually samples states from the distribution...
// Previously, I had stored each state, and thus created a distribution
// However, this is very costly, and instead, it would make more sense to store the 'observables'
// In the form of a vai
// From each step of the monte-carlo method, There is an associated observable..


//


//Obtain all observables from sample




template<class T>
class distribution
{
	T state;
public:

	distribution( const T& state ) : state( state ) {}

	// sample method is overridden ! It could call for an update... or something..
	// sample( size_t lag = 0 ) { burn( lag ); update(); return obs(  ) }
	// There will be in general some value that is a function of the state..
	// There is a function of the distribution...

	int sample( ( double )( *f_arr[] )( const T & ), double * samples, size_t len )
	{
		for( size_t i = 0; i < len; i++ )
		{
			samples[ i ] = f[ i ]( state )
		}
		return 1;
	}
};


// A distribution comes with a sample method.
// This sample method returns an array of obserables
// That are functions of the current sample
// Disfferent distributions can be sampled with different procedures
// FOr example, an MCMC is a distribution... with a sample method
// It is not an array of samples, but rather a means to obtain samples
// A distribution is a generator of samples..
// 



//each time the distribution updates, there is a 


/*********************************************

DISTRIBUTION : Object that can be sampled from

**********************************************/
/*********************************************

OBSERVABLE : Function of distribution

**********************************************/


class analyst{
private:
	double* _sample
	size_t 	n_samples;

	double 	m1;
	double 	m2;
	double 	var;

public:
	double 	compute_mean( double[], size_t );
	double 	compute_variance( double[], size_t );
	int 	analyze( double[], size_t );

};

double analyst::compute_mean( double sample[], size_t n_samples )
{
	_sample = sample;
	double acc = 0;
	for( size_t i = 0; i < n_samples; i++ )
	{
		acc += sample[ i ];
	}

	return acc / n_samples;
}


double analyst::compute_variance( double sample[], size_t n_samples )
{
	_sample = sample;
	m1 			= mean( sample, n_samples )
	double acc	= 0;
	for( size_t i = 0; i < n_samples; i++ )
	{
		acc += _sample[ i ] * _sample[ i ];
	}
	return ( acc / n_samples ) - mu * mu;
}

int analyzer::analyze( double sample[], size_t n_samples )
{
	compute_variance( sample[], n_samples );
}


int main()
{
	distribution<T> dist;
	std::vector<observable> obs;
	for( size_t i = 0; i < n_samples; i++ )
	{
		sample[ i ] = dist.sample();
	}

	analyst a;
	a.analyze( sample, n_samples );


	return 1;
}

//post-processing data for 