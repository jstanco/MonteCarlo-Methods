//created by John Stanco 6.24.18

#include "../include/help.hpp"
#include "../include/lattice.hpp"


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
	return test_lattice( 1 );
}	