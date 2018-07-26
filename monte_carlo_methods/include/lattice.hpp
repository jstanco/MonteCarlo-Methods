 /*
 * lattice.hpp

 *
 *  Created on: Apr 8, 2017
 *      Author: johnstanco
 */

#include "world"
#include <sstream>
#include </usr/local/include/armadillo>

#ifndef lattice_hpp
#define lattice_hpp


template <class T>
class tensor{

private:
  const size_t dim;
  const size_t N;
  const arma::uvec& n_max;
  std::vector<T> data;
public:
  tensor( const arma::uvec& n_max   ) : 
    dim   ( n_max.size()            ) , 
    N     ( prod( n_max )           ) ,
    n_max ( n_max                   ) ,
    data  ( std::vector<T>( prod( n_max ) ) ) {}

  T& operator()( const arma::uvec& indices )
  {
    size_t index = 0;
    size_t mult  = 1;
    for ( int i = dim - 1; i >= 0; i-- )
    {
      index  += indices( i ) * mult;
      mult   *= n_max( i );
    }
    return data[ index ];
  }

  //allows indexing of underlying array ( more efficient for iterating over all elements ).
  T& operator[](const size_t& index) { return data[ index ]; }
};


template<class T> class image;
template<class T> class lattice;

class bead
{
  friend class image<bead>;
  friend class lattice<bead>;
private:
  size_t        image_index;
  image<bead> * current_image;
  bead        * prev;
  bead        * next;
  arma::vec     pos; //in cartesian coords.

public:

  const double mass;
  const double r_c;
  const double eps;

  bead( const double mass, const double r_c, const double eps ) : 
    mass  ( mass ) ,
    r_c   ( r_c  ) ,
    eps   ( eps  ) {}

  void move( const arma::vec &r ) { pos = r; }
  arma::vec& get_pos() { return pos; }
};



template<class T>
class image
{
private:
  std::vector< std::vector< T *> >  atom_storage;
  std::vector<image *>              nn;
  size_t                            n_slices;

  //swap image to end
  int swap( const size_t im_ind, const size_t t_ind )
  {
    atom_storage[ t_ind ][ im_ind ]                   = atom_storage[ t_ind ][ atom_storage.size() - 1 ];
    atom_storage[ t_ind ][ im_ind ]->image_index      = im_ind;
    atom_storage[ t_ind ][ atom_storage.size() - 1 ]  = nullptr;

    return 1;
  }
public:
  image( const size_t n_slices = 1 ) : 
    atom_storage  ( std::vector< std::vector<T *> >() ) , 
    nn            ( std::vector<image *>()               ) , 
    n_slices      ( n_slices ) { atom_storage.resize( n_slices ); }

  void add( T* new_bead, const size_t t_ind )
  {
    atom_storage[ t_ind ].push_back( new_bead );
    new_bead->current_image = this;
    new_bead->image_index   = atom_storage[ t_ind ].size() - 1;
  }

  bead* remove( const T * old_bead, const size_t t )
  {
    /*
    if ( *old_bead->image != *this )
    {
      throw("Attempted to remove bead from wrong image!");
    }
    */
    swap( old_bead->image_index, t );
    old_bead->image        = nullptr;
    old_bead->image_index  = 0;
    return old_bead;
  }
};


std::ofstream ofs("../data/lattice.dat");


class material{
private:

  

public:
  const arma::mat&   vecs;
  const arma::mat&   atom_pos;
  const arma::vec&   masses;
  const arma::vec&   sigma;
  const arma::vec&   eps;
  const size_t       dim;
  const size_t       n_part;

  material( const arma::mat& vecs,
            const arma::mat& atom_pos,
            const arma::vec& masses,
            const arma::vec& sigma,
            const arma::vec& eps ) : 
      vecs    ( vecs        ) ,
      atom_pos( atom_pos    ) ,
      masses  ( masses      ) ,
      sigma   ( sigma       ) ,
      eps     ( eps         ) ,
      dim     ( vecs.n_cols ) ,
      n_part  ( atom_pos.n_cols ) {}

  class builder
  {
  private:
    arma::mat   _vecs;
    arma::mat   _atom_pos;
    arma::vec   _masses;
    arma::vec   _sigma;
    arma::vec   _eps;

  public:
    builder() {}
    builder& set_vecs   ( const arma::mat& vecs   )  { _vecs     = vecs;   return *this; }
    builder& set_pos    ( const arma::mat& pos    )  { _atom_pos = pos;    return *this; }
    builder& set_masses ( const arma::vec& masses )  { _masses   = masses; return *this; }
    builder& set_sigma  ( const arma::vec& sigma  )  { _sigma    = sigma;  return *this; }
    builder& set_eps    ( const arma::vec& eps    )  { _eps      = eps;    return *this; }
    material build() { return material( _vecs, _atom_pos, _masses, _sigma, _eps ); }
  };
};


template<class T>
class lattice
{
private:
  const size_t       dim;
  const size_t       n_images;
  const arma::uvec & r_max;  //in lattice coords ( n_1, ..., n_d ), n_i an integer. 
  const material   & matl;
  std::vector<T *>   beads;  //redundant pointers, use for selecting beads uniformly.
  tensor<image<T> *> images; //image should not just store beads, but store each time slice separately...

  T* init_bead( const size_t index, const arma::uvec& indices, const size_t p_index, const size_t t )
  {

    T * new_bead  = new bead( matl.masses( p_index ), matl.sigma( p_index ), matl.eps( p_index ) );
    arma::vec r   = matl.vecs * indices + matl.atom_pos.col( p_index );
    new_bead->move( r );
    images[ index ]->add( new_bead, t );
    beads.push_back( new_bead );
    return new_bead;

  }


  int fill_nn( const size_t index, const arma::uvec& indices )
  {
    image<T> * current_image = images[ index ];

    size_t mult = 1;
    for ( size_t i = 0; i < dim; i++ )
    { 
      size_t j                    = indices( i );
      size_t dn                   = mult * ( ( j + 1 ) % r_max( i ) - j );
      current_image->nn[ 2 * i ]  = images[ index + dn ];
      mult *= r_max( i );
    }
    return 1;
  }


  arma::uvec indices( const size_t index )
  {
      arma::uvec ind( dim );
      size_t loop_max = dim - 1;
      size_t y_j  = 0;
      size_t mult = 1;
      for ( size_t j = 0; j < loop_max; j++ )
      {
        ind( loop_max - j )  = ( ( index - y_j ) / mult ) % r_max( loop_max - j );
        y_j                 += ind( loop_max - j ) * mult;
        mult                *= r_max( loop_max - j );

      }
      ind( 0 ) = ( ( index - y_j ) / mult );
      return ind;
  }
  

  int populate_image( const size_t index, const size_t n_slices )
  {

    images[ index ] = new image<T>( n_slices );
    arma::uvec ind  = indices( index );
    for ( size_t i = 0; i < matl.n_part; i++ )
    {
      
      T * head = init_bead( index, ind, i, 0 );
      T * tmp  = head;

      for ( size_t j = 1; j < n_slices; j++ )
      {
      
        T * new_bead    = init_bead( index, ind, i, j );   //have to take into account time slices

        new_bead->prev  = tmp;                               
        tmp->next       = new_bead;
        tmp             = new_bead;

      }

      tmp->next = head;
      head->prev = tmp;

      for ( size_t j = 0; j < dim; j++ )
      {
        ofs << head->get_pos()( j ) << "\t";
      }
      ofs << std::endl;
    }

    return 1;
  }


  int fill( const size_t n_slices = 1 )
  {
    for ( size_t i = 0; i < n_images; i++ )
    {
      populate_image( i, n_slices );
    }
    return 1;
  }

  
public:
  
  /*
  int move_iter( image<T> * iter, const arma::uvec & dn )
  {
    //iterate through lattice
    for ( size_t i = 0; i < dim; i++ )
    {
      size_t move_dir = 0;
      
      if ( dn( i ) < 0 )
      {
        move_dir = 1;
      }
      
      for ( size_t j = 0; j < dn( i ); j++ )
      {
        size_t nn_index = i + move_dir;
        if( iter->nn[ nn_index ] )
        {
          iter = iter->nn[ nn_index ];
        }
        else 
        {
          break;
        }
      }
    }

    return 0;
  }
  */

  int delete_beads()
  {
    for ( size_t i = 0; i < beads.size(); i++ )
    {
      if( beads[ i ] ) delete beads[ i ];
    }
    return 1;
  }


  int delete_images()
  {
    for ( size_t i = 0; i < n_images; i++ )
    {
      if( images[ i ] ) delete images[ i ];
    }
    return 1;
  }
  
  
  int update_image( T * bead, const arma::vec & r )
  {
    arma::uvec im_index  = floor( inv( matl.vecs ) * r );  //consider storing inv( vecs ) if too costly  
    bead->current_image  = images( im_index );
    return bead->move( r );
  }


  lattice( const material& matl, const arma::uvec& r_max, const size_t n_slices = 1 ) :
    matl    ( matl            ) ,
    dim     ( matl.dim        ) , 
    r_max   ( r_max           ) , 
    n_images( prod( r_max )   ) ,
    images  ( tensor<image<T> *>( r_max ) ) { fill( n_slices ); }
  

  ~lattice() { delete_beads(); delete_images(); }

};


#endif /* lattice_hpp */
