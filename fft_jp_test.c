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
#define N 64
#define M 32
#define L 2*M_PI
#define W_REAL 0

int main( int argc, char **argv )
{
  vector *in, *out, *grid;
  FILE *fp = fopen("fft_jp.dat", "w");
  if( !fp )
    EXIT_WITH_PERROR("file open failed in main: ")

  create_vector( &in, N );
  create_vector( &out, N );
  create_vector( &grid, N );

  vec_fill_grid( grid, L );
  vec_fill_cosine( in, M, L );

  fft( in, out, N );

  write_data( fp, out, grid, N, W_REAL );

  free( in );
  free( out );
  free( grid );

  fclose( fp );


  return 0;
}
