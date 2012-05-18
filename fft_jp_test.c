/*
 * fft_jp_test.c
 *
 * Main test program for my fft_jp code.
 *
 * Author: Jeffrey Picard
 */
/* Standard Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* My headers */
#include "fft_jp.h"

/* Constants */
#define N 524288
#define M N
#define L 2*M_PI

int main( int argc, char **argv )
{
  vector *in, *out, *grid;
  const int n = N;
  const int m = M;
  FILE *fp = fopen("fft_jp.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  create_vector( &in, n );
  create_vector( &out, n );
  create_vector( &grid, n );

  vec_fill_grid( grid, L );
  vec_fill_cosine( in, n, L );

  fft( in, out, n );

  write_data( fp, out, grid, n, W_REAL );

  free( in );
  free( out );
  free( grid );

  fclose( fp );

  return 0;
}
