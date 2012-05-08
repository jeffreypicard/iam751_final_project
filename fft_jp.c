/*
 * fft_jp.c
 *
 * A parallel fft implementation written as a final project
 * for IAM751 Spring 2012.
 *
 * Author: Jeffrey Picard (jpicardnh@gmail.com)
 */

/* Standard headers */
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>

/* My headers */
#include "fft_jp.h"

/*
 * fft
 *
 * Takes an input vector, output vector, number of points
 * and w^1 = -i.
 * Performs a Fast Foruier Transform on the input vector
 * and stores it in the output vector.
 */
void fft( vector *in, vector *out, int n, complex double w )
{
  int k;
  if( n == 1 )
    out[0] = in[0];
  else
  {
    vector *x, *y, *p, *q;
    create_vector( &x, n/2 );
    create_vector( &y, n/2 );
    create_vector( &p, n/2 );
    create_vector( &q, n/2 );
    for( k = 0; k < n/2; k++ )
    {
      //x[k] = in[2*k]
      //y[k] = in[2*k + 1]
      VEC( x, k ) = VEC( in, 2*k );
      VEC( y, k ) = VEC( in, 2*k + 1 );
    }
    fft( x, p, n/2, w*w );
    fft( y, q, n/2, w*w );
    for( k = 0; k < n; k++ )
    {
      //out[k] = p[ k % (n/2) ] + pow(w,k)*q[k % (n/2) ];
      VEC( out, k ) = VEC( p, k % (n/2) ) + pow(w,k)*VEC( q, k % (n/2) );
    }
    destroy_vector( x );
    destroy_vector( y );
    destroy_vector( p );
    destroy_vector( q );
  }
}

/*
 * fill_grid_vec
 *
 * Takes a vector and a scalar and fills the vector
 * to be an evenly spaced grid from 0 ... scale.
 */
void fill_grid_vec( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = i * scale / v->n;
}

/*
 * create_vector
 *
 * Takes a pointer to a pointer to a vector and a size.
 * Allocates space for the vector struct and the number
 * of specified values.
 */
void create_vector( vector **v, int n )
{
  *v = malloc( sizeof(vector) );
  if( !*v )
    EXIT_WITH_PERROR("malloc failed in create_vector: ")
  (*v)->vals = calloc( sizeof(val_type), n );
  if( !(*v)->vals )
    EXIT_WITH_PERROR("malloc failed in create_vector: ")
  (*v)->n = n;
}

/*
 * destroy_vector
 *
 * Takes a vector pointer and frees all of its
 * associated memory, including the struct.
 */
void destroy_vector( vector *v )
{
  free( v->vals );
  free( v );
}
