/*
 * fftw_test.c
 *
 * Testing code using the fftw library./
 *
 * Author: Jeffrey Picard
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>

#include <fftw.h>

/* Constants */
#define N 1000

int main( int argc, char **argv )
{
  fftw_complex in[N], out[N];
  fftw_plan p;

  p = fftw_create_plan( N, FFTW_FORWARD, FFTW_ESTIMATE );

  fftw_one( p, in, out );

  fftw_destroy_plan( p );

  return 0;
}
