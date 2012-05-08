/*
 * fftw_test.c
 *
 * Testing code using the fftw library.
 *
 * Author: Jeffrey Picard (jpicardnh@gmail.com)
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 1000
#define L 2*M_PI
#define W_REAL 0

/* Function Prototypes */
void fftw_complex_fill_sine( fftw_complex *, const int );
void fftw_complex_fill_cosine( fftw_complex *, const int );
void fill_grid( double *, const int, double );
void fftw_write_data( FILE*, fftw_complex *, double *, const int, int );

int main( int argc, char **argv )
{
  fftw_complex in[N], out[N];
  double grid[N];
  fftw_plan p;
  FILE *fp = fopen("fftw.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  fill_grid( grid, N, L );
  fftw_complex_fill_cosine( in, N );

  p = fftw_create_plan( N, FFTW_FORWARD, FFTW_ESTIMATE );

  fftw_one( p, in, out );

  fftw_write_data( fp, out, grid, N, W_REAL );

  fftw_destroy_plan( p );

  fclose( fp );

  return 0;
}

/*
 * fftw_write_data
 *
 * Takes a FIEL*, an fftw_complex vector a double* grid vector, 
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
    fprintf( fp, "%g %g\n", grid[i], v[i].re );
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
    v[i].re = sin( i * L / n );
    v[i].im = 0.;
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
    v[i].re = cos( i * L / n );
    v[i].im = 0.;
  }
}
