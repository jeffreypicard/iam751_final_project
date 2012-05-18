/*
 * fftw_timing.c
 *
 * Timing code using the fftw library.
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
#define N 524288
#define M N
#define L 2*M_PI
#define W_REAL 0
#define MAX_TESTS 100

/* Function Prototypes */
void fftw_complex_fill_sine( fftw_complex *, const int );
void fftw_complex_fill_cosine( fftw_complex *, const int );
void fill_grid( double *, const int, double );
void fftw_write_data( FILE*, fftw_complex *, double *, const int, int );

int main( int argc, char **argv )
{
  fftw_complex *in, *out;
  fftw_plan p;
  int num_tests = 0;
  int n = 2;
  /* Just in case, might be nice to have these */
  int points_arr[MAX_TESTS];
  double time_data[MAX_TESTS];
  double t_b, t_e;
  int i;
  if( argc == 2 )
  {
    num_tests = atoi( argv[1] );
    if( num_tests > MAX_TESTS )
      EXIT_WITH_ERROR("Error: max number of tests is %d\n", MAX_TESTS );
  }
  else
    EXIT_WITH_ERROR("Usage: fftw_test <num_pows_of_two>\n");

  FILE *fp = fopen("fftw_timing.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  for( i = 0; i < num_tests; i++ )
  {
    t_b = WTime();
    in = fftw_malloc( sizeof(fftw_complex) * n );
    out = fftw_malloc( sizeof(fftw_complex) * n );

    if( !in || !out )
      EXIT_WITH_PERROR("malloc failed in main")


    fftw_complex_fill_cosine( in, n );

    //p = fftw_create_plan( N, FFTW_FORWARD, FFTW_ESTIMATE );
    p = fftw_plan_dft_1d( n, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

    /*fftw_one( p, in, out );*/
    fftw_execute( p );

    fftw_destroy_plan( p );
    fftw_free( in );
    fftw_free( out );

    t_e = WTime();

    points_arr[i] = n;
    time_data[i] = t_e - t_b;

    /*fftw_write_data( fp, out, grid, N, W_REAL );*/
    fprintf( fp, "%d %g\n", n, t_e - t_b );
    fprintf( stderr, "%d n = %d\n", i, n );

    n *= 2;
  }

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
