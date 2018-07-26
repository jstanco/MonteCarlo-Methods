//created by John Stanco 6.24.18

#include "../include/help.hpp"
#include "../include/lattice.hpp"
#include "../../Signal/include/timer.hpp"


int test_2d_loop( size_t nx, size_t ny )
{
	for ( size_t i = 0; i < nx; i++ )
	{
		for ( size_t j = 0; j < ny; j++ )
		{
			size_t k = i * ny + j;
			std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << std::endl;
		}
	}

	std::cout << std::endl;

	for( size_t k = 0; k < nx * ny; k++ )
	{

		size_t j = k % ny;
		size_t i = ( k - j ) / ny;

		std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << std::endl;
	}
	return 1;
}


int test_3d_loop( size_t nx, size_t ny, size_t nz )
{
	for ( size_t i = 0; i < nx; i++ )
	{
		for ( size_t j = 0; j < ny; j++ )
		{
			for ( size_t k = 0; k < nz; k++ )
			{
				size_t l = i * ny * nz + j * nz + k;
				std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tl: " << l << std::endl;
			}
			
		}
	}

	std::cout << std::endl;

	for( size_t l = 0; l < nx * ny * nz; l++ )
	{

		size_t k = l % nz;
		size_t j = ( ( l - k) / nz) % ny;
		size_t i = ( l - j * nz - k ) / ( ny * nz );

		std::cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tl: " << l << std::endl;
	}
	return 1;
}


template <class T>
T product( const arma::Col<T>& x )
{
	T prod = 1; 
	for( size_t i = 0; i < x.size(); i++ )
	{
		prod *= x( i );
	}
	return prod;
}


arma::uvec compute_stride( const arma::uvec& r_max )
{
  size_t dim = r_max.n_elem;
  arma::uvec stride_len( dim );
  stride_len( dim - 1 ) = 1;
  for( int i = dim - 2; i >= 0; i-- )
  {
    stride_len( i ) = stride_len( i + 1 ) * r_max( i + 1 );
  }
  return stride_len;
}


int to_index( const arma::uvec& indices, const arma::uvec& stride_len )
{
	size_t dim = indices.size();
	size_t ind = 0;
	for ( int i = dim - 1; i >= 0; i-- )
	{
		ind  += indices( i ) * stride_len( i );
	}	
  return ind;
}




int to_index_vec_mult( const arma::uvec& indices, const arma::uvec& stride_len )
{
  return dot( stride_len, indices );
}


arma::uvec to_vec_indices( const arma::uvec& ns, const size_t index )
{
    size_t dim      = ns.n_elem;
    size_t y_i      = 0;
    size_t mult     = 1;
    arma::uvec xs( dim );
    for ( int i = dim - 1; i > 0; i-- )
    {

      xs( i )  = ( ( index - y_i ) / mult ) % ns( i );
      y_i     += xs( i ) * mult;
      mult    *= ns( i );

    }

    xs( 0 ) = ( ( index - y_i ) / mult );
    return xs;
}


//testing indexing over n - dimensional loop
int test_nd_loop( const arma::uvec& ns )
{

	size_t dim 		= ns.size();
	size_t n_tot	= prod( ns );
	size_t loop_max = dim - 1;

	for ( size_t l = 0; l < n_tot; l++ )
	{
    arma::uvec xs = to_vec_indices( ns, l );

		for ( size_t i = 0; i < dim; i++ )
		{
			std::cout << i << ": " << xs( i ) << "\t";
		}

		std::cout << "l: " << l;
		std::cout << "\tl: " << to_index( xs, compute_stride( ns ) ) << std::endl;
	}

	return 1;
}


//random integer indices
arma::uvec rand_indices( const arma::uvec & r_max )
{
  arma::uvec indices( r_max.n_elem );
  for ( size_t i = 0; i < r_max.n_elem; i++ )
  {
    indices( i ) = rand_int( r_max( i ) );
  }
  return indices;
}


int test_index_time( const arma::uvec& r_max )
{
  size_t i1 = 0;
  size_t i2 = 0;
  timer t1, t2;
  arma::uvec stride_len = compute_stride( r_max );
  for( size_t i = 0; i < 100000; i++ )
  {
    arma::uvec indices = rand_indices( r_max );
    t1.start();
    i1 = to_index( indices, stride_len );
    t1.stop();

    t2.start();
    i2 = to_index_vec_mult( indices, stride_len );
    t2.stop();

    if( i1 != i2 ) return 0;
  }

  std::cout << t1.t_sec() / t2.t_sec() << std::endl;

	return 1;
}


int test_loop(){
	size_t d = 2; //d-dimensional lattice
	arma::uvec ns( d );
	for ( size_t i = 0; i < d; i++ )
	{
		ns( i ) = 2;
	}

	return test_nd_loop( ns );
}


arma::mat bravais_lat_3d( int i_brav, const arma::vec& spacing, const arma::vec& angles )
{
  arma::mat vecs( 3, 3 );

  double c4;
  double tx;
  double ty;
  double tz;
  double u;
  double v;
  double v2y;
  double v2z;

  double a = spacing( 0 );
  double b = spacing( 1 );
  double c = spacing( 2 );

  double alpha = angles( 0 );
  double beta  = angles( 1 );
  double gamma = angles( 2 );

  switch(i_brav)
  {
    case 1 :                            //cubic P (sc)
      vecs  << 1 << 0 << 0 << arma::endr
            << 0 << 1 << 0 << arma::endr
            << 0 << 0 << 1 << arma::endr;
      vecs *= a;
      break;
    case 2 :                            //cubic F (fcc)
      vecs  << -1 << 0 << 1 << arma::endr
            << 0  << 1 << 1 << arma::endr
            << -1 << 1 << 0 << arma::endr;
      vecs *= ( a / 2 );
      break;
    case 3 :                            //cubic I (bcc)
      vecs  << 1  << 1  << 1 << arma::endr
            << -1 << 1  << 1 << arma::endr
            << -1 << -1 << 1 << arma::endr;
      vecs *= ( a / 2 );
      break;
    case 4 :                            //Hex & Trigonal P
        vecs  << 1   << 0             << 0     << arma::endr 
        	  << -.5 << sqrt( 3 ) / 2 << 0     << arma::endr 
        	  << 0   << 0             << c / a << arma::endr;
        vecs *= a;
      break;
    case 5 :                            //Trigonal R (3 - fold axis around <111>)
      c4 = cos( alpha );
      tx = sqrt( ( 1 - c4 ) / 2 );
      ty = sqrt( ( 1 - c4 ) / 6 );
      tz = sqrt( ( 1 + 2 * c4 ) / 3 );
      vecs  << tx << -ty << tz << arma::endr
            << 0 << 2 * ty << tz << arma::endr
            << -tx << -ty << tz << arma::endr;
      vecs *= a;
      break;
    case -5 :                           //Trigonal R (3 - fold axis around z)
      c4 = cos( alpha );
      tx = sqrt( ( 1 - c4 ) / 2 );
      ty = sqrt( ( 1 - c4 ) / 6 );
      tz = sqrt( ( 1 + 2 * c4 ) / 3 );
      u = tz - 2 * sqrt( 2 ) * ty;
      v = tz + sqrt( 2 ) * ty;
      vecs  << u << v << v << arma::endr
            << v << u << v << arma::endr
            << v << v << u << arma::endr;
      vecs *= ( a / sqrt( 3 ) );
      break;
    case 6 :                            //Tetragonal P (st)
      vecs  << 1 << 0 << 0 << arma::endr
            << 0 << 1 << 0 << arma::endr
            << 0 << 0 << c / a << arma::endr;
      vecs *= a;
      break;
    case 7 :                            //Tetragonal I (bct)
      vecs  << 1 << -1 << c / a << arma::endr
            << 1 << 1 << c / a << arma::endr
            << -1 << -1 << c / a << arma::endr;
      vecs *= (a / 2);
      break;
    case 8 :                            //Orthorhombic P
      vecs  << a << 0 << 0 << arma::endr
            << 0 << b << 0 << arma::endr
            << 0 << 0 << c << arma::endr;
      break;
    case 9 :                            //Orthorhombic base-centered (bco)
      vecs  << a / 2 << b / 2 << 0 << arma::endr
            << -a / 2 << b / 2 << 0 << arma::endr
            << 0 << 0 << c << arma::endr;
      break;
    case -9 :                           //Orthorhombic base-centered (alternate descr.)
      vecs  << a / 2 << -b / 2 << 0 << arma::endr
            << a / 2 << b / 2 << 0 << arma::endr
            << 0 << 0 << c << arma::endr;
      break;
    case 10 :                           //Orthorhombic face-centered
      vecs  << a / 2 << 0 << c / 2 << arma::endr
            << a / 2 << b / 2 << 0 << arma::endr
            << 0 << b / 2 << c / 2 << arma::endr;
      break;
    case 11 :                           //Orthorhombic body-centered
      vecs  << a / 2 << b / 2 << c / 2 << arma::endr
            << -a / 2 << b / 2 << c / 2 << arma::endr
            << -a / 2 << -b / 2 << c / 2 << arma::endr;
      break;
    case 12 :                           //Monoclinic P unique axis c
      vecs  << a << 0 << 0 << arma::endr
            << b * sin( gamma ) << b * cos( gamma ) << 0 << arma::endr
            << 0 << 0 << c << arma::endr;
      break;
    case -12 :                          //Monoclinic P unique axis b
      vecs  << a << 0 << 0 << arma::endr
            << 0 << b << 0 << arma::endr
            << c * cos( beta ) << 0 << c * sin( beta ) << arma::endr;
      break;
    case 13 :                           //Monoclinic base-centered
      vecs  << a / 2 << 0 << -c / 2 << arma::endr
            << b * cos( gamma ) << b * sin( gamma ) << 0 << arma::endr
            << a / 2 << 0 << c / 2 << arma::endr;
      break;
    case 14 :                           //Triclinic
      v2y = c * ( cos( alpha ) - cos( beta ) * cos( gamma ) ) / sin( alpha );
      v2z = c * sqrt( 1 + 2 * cos( alpha ) * cos( beta ) * cos( gamma )
                     - pow( cos( alpha ), 2 ) - pow( cos( beta ), 2 ) - pow( cos( gamma ), 2 ) )
      / sin( gamma );
      vecs  << a << 0 << 0 << arma::endr
            << b * cos( gamma ) << b * sin( gamma ) << 0 << arma::endr
            << c * cos( beta ) << v2y << v2z << arma::endr;
      break;
    default : std::cout << i_brav << " is not a valid Bravais lattice index" << std::endl;
      break;
  }
  return trans( vecs );
}


arma::mat fcc_coords( const arma::mat& coords, const double a )
{
	//converts a primitive basis into fcc coordinates
  return arma::mat();
}


//gen_sym

arma::mat gen_sym(  );



double toRad( double theta_deg )
{
	return theta_deg * pi / 180;
}


int test_lattice( int i_brav )
{
	//create a lattice

	clock_t t = clock();

	size_t      dim = 3;
	size_t 		n_part = 8;

	double 		a = 1e-10;
	double		b = 1e-10;
	double 		c = 2e-10;

	double alpha = pi / 3;
	double beta  = 2 * pi / 3;
	double gamma = pi / 2;

	arma::vec spacing;
	arma::vec angles;

	spacing << a << b << c << arma::endr;
	angles  << alpha << beta << gamma << arma::endr;

  	arma::mat  vecs 		= bravais_lat_3d( i_brav, spacing, angles );
  	arma::mat  atom_pos;

  	atom_pos 	<< 0 << 2 << 2 << 0 << 3 << 3 << 1 << 1 << arma::endr
  				<< 0 << 2 << 0 << 2 << 3 << 1 << 3 << 1 << arma::endr
  				<< 0 << 0 << 2 << 2 << 3 << 1 << 1 << 3 << arma::endr;

  	atom_pos *= a / 4;



  	arma::uvec r_max 		= arma::ones<arma::uvec>( dim ) * 20;
  	arma::vec  masses 	= arma::ones( n_part );
  	arma::vec  sigma 		= arma::ones( n_part ) * a;
  	arma::vec  eps 			= arma::ones( n_part );


  	/*material matl = material::builder().set_vecs	( vecs )
										.set_pos	( atom_pos )
										.set_masses	( masses )
										.set_sigma	( sigma )
										.set_eps	( eps )
										.build();*/



  	material 		matl( vecs, atom_pos, masses, sigma, eps );
  	t = clock() - t;

  	std::cout << "\nTime to specify material:\t" << ( double )t / CLOCKS_PER_SEC << std::endl;


  	t = clock();
  	lattice<bead> 	lat( matl, r_max, 1 );
  	t = clock() - t;

  	std::cout << "\nTime to construct lattice:\t" << ( double )t / CLOCKS_PER_SEC << "\n" << std::endl;
  	return 1;
}


int main()
{
	//return test_lattice( 1 );
  arma::uvec ns( 3 );
  ns << 20 << 17 << 5 << arma::endr;
  std::cout << test_index_time( ns ) << std::endl;
  return 1;
}	