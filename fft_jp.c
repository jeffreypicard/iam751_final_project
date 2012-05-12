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
#include <mpi.h>

/* My headers */
#include "fft_jp.h"

/*
 * fft
 *
 * Takes an input vector, output vector and a number of points
 * Performs a Fast Fourier Transform on the input vector
 * and stores it in the output vector.
 */
void fft( vector *in, vector *out, int n )
{
  int k;
  if( n == 1 )
    VEC( out, 0 ) = VEC( in, 0 );
  else
  {
    complex double wn = cexp( -2*M_PI*I / n );
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
    fft( x, p, n/2 );
    fft( y, q, n/2 );
    for( k = 0; k < n; k++ )
    {
      //out[k] = p[ k % (n/2) ] + pow(w,k)*q[k % (n/2) ];
      VEC( out, k ) = VEC( p, k % (n/2) ) + cpow(wn,k)*VEC( q, k % (n/2) );
    }
    destroy_vector( x );
    destroy_vector( y );
    destroy_vector( p );
    destroy_vector( q );
  }
}

/*
 * fft_mpi
 *
 * Takes an input vector, output vector and a number of points
 * Performs a Fast Fourier Transform on the input vector
 * and stores it in the output vector.
 * 
 * Uses MPI for parallelization 
 */
void fft_mpi( vector *in, vector *out, int n )
{
  int k;
  if( n == 1 )
    VEC( out, 0 ) = VEC( in, 0 );
  else
  {
    complex double wn = cexp( -2*M_PI*I / n );
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
    fft( x, p, n/2 );
    fft( y, q, n/2 );
    for( k = 0; k < n; k++ )
    {
      //out[k] = p[ k % (n/2) ] + pow(w,k)*q[k % (n/2) ];
      VEC( out, k ) = VEC( p, k % (n/2) ) + cpow(wn,k)*VEC( q, k % (n/2) );
    }
    destroy_vector( x );
    destroy_vector( y );
    destroy_vector( p );
    destroy_vector( q );
  }
}

/*
 * vec_fill_sine
 *
 * Takes a vector and a scalar and fills the vector
 * with sine values.
 */
void vec_fill_sine( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = sin( i * scale / v->n );
}

/*
 * vec_fill_cosine
 *
 * Takes a vector and a scalar and fills the vector
 * with cosine values.
 */
void vec_fill_cosine( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = cos( i * scale / v->n );
}

/*
 * vec_fill_grid
 *
 * Takes a vector and a scalar and fills the vector
 * to be an evenly spaced grid from 0 ... scale.
 */
void vec_fill_grid( vector *v, val_type scale )
{
  int i;
  for( i = 0; i < v->n; i++ )
    VEC( v, i ) = (double)i / (double)v->n;
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
  (*v) = malloc( sizeof(vector) );
  if( !*v )
    EXIT_WITH_PERROR("malloc failed in create_vector: ")
  (*v)->vals = calloc( n, sizeof(val_type) );
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
  if( !v ) return;
  free( v->vals );
  free( v );
}

/*
 * write_data
 *
 * Takes a FILE*, a data vector, a grid vector 
 * a size and a mode.
 * mode 0: write real values
 * mode 1: write complex values
 * mode 2: write both real and comples values
 */
void write_data( FILE *fp, vector *v, vector *grid, 
                      const int n, int mode )
{
  int i;
  for( i = 0; i < n; i++ )
    fprintf( fp, "%g %g\n", creal(VEC( grid, i )), creal(VEC( v, i )) );
}
