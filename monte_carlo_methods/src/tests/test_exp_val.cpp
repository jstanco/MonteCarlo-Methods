//created by John Stanco 7.1.2018


#include "../include/help.hpp"


double e_val_1( double ( * f )( double ), double dist[], size_t len )
{
	double e_val = 0;
	for( size_t i = 0; i < len; i++ )
	{
		e_val = ( e_val * i + f( dist[ i ] ) ) / ( i + 1 );
	}
	return e_val;
}


double e_val_2( double ( * f )( double ), double dist[], size_t len )
{
	double sum = 0;
	for( size_t i = 0; i < len; i++ )
	{
		sum += f( dist[ i ] );
	}
	return sum / len;
}


int build_dist( double dist[], size_t len )
{	
	for ( size_t i = 0; i < len; i++ )
	{
		dist[ i ] = pow( exp( - i ), .1 );
	}
	return 1;
}


double eye( double x )
{
	return x;
}


double test_f( double x )
{
	return 3 * x + 5;
}


int test_exp_val()
{
	size_t len = 25;
	double dist[ len ];

	build_dist( dist, len );

	double exp_1 = e_val_1( test_f, dist, len );
	double exp_2 = e_val_2( test_f, dist, len) ;
	return ( exp_1 == exp_2 );
}

//compute observables based on array of function pointers -> create framework that allows for classes



int main()
{
	std::cout << test_exp_val() << std::endl;
}
