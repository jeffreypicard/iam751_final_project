/*
 * fftw_2d_test.c
 *
 * Testing code using the fftw library for the
 * 2d fft.
 *
 * Author: Jeffrey Picard (jpicardnh@gmail.com)
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 1024
#define M N*N
#define L 2*M_PI
#define W_REAL 0

/* Function Prototypes */
void fftw_complex_fill_sine( fftw_complex *, const int );
void fftw_complex_fill_cosine( fftw_complex *, const int );
void fftw_complex_fill_cosine_2d( fftw_complex *, const int, const int);
void fill_grid( double *, const int, double );
void fftw_write_data( FILE*, fftw_complex *, double *, const int, int );

int main( int argc, char **argv )
{
  fftw_complex *in, *out;
  double *grid;
  fftw_plan p;
  FILE *fp = fopen("fftw_2d.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  grid = calloc( sizeof(double), M );
  in = fftw_malloc( sizeof(fftw_complex) * M );
  out = fftw_malloc( sizeof(fftw_complex) * M );

  if( !grid || !in || !out )
    EXIT_WITH_PERROR("malloc failed in main")


  fill_grid( grid, M, L );
  //fftw_complex_fill_cosine( in, M );
  fftw_complex_fill_cosine_2d( in, N, N );
  /*int i;
  for( i = M; i < N; i++ )
    in[i] = 0;*/

  //p = fftw_create_plan( N, FFTW_FORWARD, FFTW_ESTIMATE );
  //p = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  p = fftw_plan_dft_2d( N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

  //fftw_one( p, in, out );
  fftw_execute( p );

  fftw_write_data( fp, out, grid, M, W_REAL );

  fftw_destroy_plan( p );
  fftw_free( in );
  fftw_free( out );
  free( grid );

  fclose( fp );

  return 0;
}

/*
 * fftw_write_data
 *
 * Takes a FILE*, an fftw_complex vector a double* grid vector, 
 * a size and a mode.
 * mode 0: write real values
 * mode 1: write complex values
 * mode 2: write both real and comples values
 */
void fftw_write_data( FILE *fp, fftw_complex *v, double *grid, 
                      const int n, int mode )
{
  int i;
  for( i = 0; i < n; i++ )
    fprintf( fp, "%g %g\n", grid[i], creal(v[i]) );
}

/*
 * fill_grid
 *
 * Takes a double vector a size and a scalar and fills
 * the vector evenly from 0 ... scale
 */
void fill_grid( double *v, const int n, double scale )
{
  int i;
  for( i = 0; i < n; i++ )
    v[i] = (double)i / (double)n;
}

/*
 * fftw_complex_fill_sine
 *
 * Takes an fftw_complex vector and the size of the vector
 * and fills it with sine values.
 */
void fftw_complex_fill_sine( fftw_complex *v, const int n )
{
  int i;
  for( i = 0; i < n; i++ )
  {
    /*v[i].re = sin( i * L / n );
    v[i].im = 0.;*/
    v[i] = sin( i * L / n );
  }
}

/*
 * fftw_complex_fill_cosine
 *
 * Takes an fftw_complex vector and the size of the vector
 * and fills it with cosine values.
 */
void fftw_complex_fill_cosine( fftw_complex *v, const int n )
{
  int i;
  for( i = 0; i < n; i++ )
  {
    /*v[i].re = cos( i * L / n );
    v[i].im = 0.;*/
    v[i] = cos( i * L / n );
    /*v[i] = cos( i * L / 10 );*/
  }
}

/*
 * fftw_complex_fill_cosine_2d
 *
 * Takes an fftw_complex vector and the size of the vector
 * and fills it with cosine values in two dimensions using 
 * row major order.
 */
void fftw_complex_fill_cosine_2d( fftw_complex *v, const int n, const int m )
{
  int i, j;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < m; j++ )
      v[(i*m) + j] = cos( j * L / m );
  }
}
